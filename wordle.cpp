#include <fstream>
#include <iostream>
#include <math.h>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// A set of words from a WordList, represented as a sorted list of indices that
// are in the set.
using WordSet = std::vector<int>;

// A list of words, using only characters a-z, all of the same length.
struct WordList {
  std::vector<std::string> words;
  int num_letters = 0;

  // Return a set representing all the words in this list.
  WordSet as_set() const {
    WordSet set;
    set.reserve(words.size());
    for (int i = 0; i < words.size(); i++) {
      set.push_back(i);
    }
    return set;
  }
};

// Possible per-letter guess outcomes.
enum Outcome {
  // Right letter in right location.
  EXACT_MATCH,
  // Right letter in wrong location.
  PARTIAL_MATCH,
  // Wrong letter.
  NO_MATCH,
};

// Guess outcomes for a full word. Response[i] is the outcome for guess[i].
using Response = std::vector<Outcome>;

// Print the given error message and abort.
void die(const std::string& message) {
  std::cout << message << std::endl;
  exit(1);
}

// Score guess against target.
Response ScoreGuess(const std::string& guess, const std::string& target) {
  if (guess.size() != target.size()) {
    die("Can't score guess and target of different sizes. guess=" + guess + ", target=" + target);
  }

  // Get all the characters in target that aren't an exact match.
  std::unordered_map<char, int> potential_partials;
  for (int i = 0; i < guess.size(); i++) {
    if (guess[i] != target[i]) {
      ++potential_partials[target[i]];
    }
  }

  // Now build the response.
  Response response(guess.size(), NO_MATCH);
  for (int i = 0; i < guess.size(); i++) {
    if (guess[i] == target[i]) {
      response[i] = EXACT_MATCH;
    } else if (potential_partials[guess[i]] > 0) {
      response[i] = PARTIAL_MATCH;
      --potential_partials[guess[i]];
    }
  }

  return response;
}

// Converts the given response into a unique integer code.
int ResponseToCode(const Response& response) {
  int code = 0;
  for (const Outcome& outcome : response) {
    code *= 3;
    switch (outcome) {
    case EXACT_MATCH:
      break;
    case PARTIAL_MATCH:
      code += 1;
      break;
    case NO_MATCH:
      code += 2;
      break;
    }
  }
  return code;
}

// Filters the given input_set (from word_list) to only those where the given
// guess would have elicited the given response.
WordSet FilterWordSet(const WordList& word_list, const WordSet& input_set,
		      const std::string& guess, const Response& response) {
  WordSet output_set;
  for (const int i : input_set) {
    const std::string& word = word_list.words[i];
    if (ScoreGuess(guess, word) == response) {
      output_set.push_back(i);
    }
  }
  return output_set;
}

double ComputeEntropy(const std::unordered_map<int, int>& distribution) {
  // First count the total number of entries in the distribution.
  int N = 0;
  for (const auto& entry : distribution) {
    N += entry.second;
  }

  // Then tally up entropy.
  double entropy;
  for (const auto& entry : distribution) {
    const double P = 1.0 * entry.second / N;
    entropy -= P * log(P);
  }

  return entropy;
}

// Display the response for this guess, using ANSI color codes.
void DisplayResponse(const std::string& guess, const Response& response) {
  if (guess.size() != response.size()) {
    die("Can't display guess and target of different sizes. guess=" + guess);
  }

  for (int i = 0; i < guess.size(); i++) {
    switch (response[i]) {
      case NO_MATCH:
	std::cout << "\033[37;40m" << guess[i] << "\033[0m";
	break;
      case EXACT_MATCH:
	std::cout << "\033[37;42m" << guess[i] << "\033[0m";
	break;
      case PARTIAL_MATCH:
	std::cout << "\033[37;43m" << guess[i] << "\033[0m";
	break;
    }
  }
  std::cout << std::endl;
}

// Print all words in the given set.
std::string WordSetToString(const WordList& word_list, const WordSet& set) {
  std::ostringstream ss;
  bool first = true;
  ss << "{";
  for (const int i : set) {
    if (!first) {
      ss << ", ";
    }
    ss << word_list.words[i];
    first = false;
  }
  ss << "}";
  return ss.str();
}

// Read a world list from the given filename.
WordList ReadWordList(const std::string& filename) {
  // Open the file.
  std::ifstream file(filename);
  if (file.fail()) {
    die("Could not open wordlist: " + filename);
  }

  // Read in the words.
  WordList list;
  std::string word;
  while (file >> word) {
    const int num_letters = word.size();
    if (list.num_letters > 0 && list.num_letters != word.size()) {
      die("Got word with wrong number of letters: " + word);
    }
    for (char c : word) {
      if (c < 'a' || c > 'z') {
	die("Got invalid character in word: " + word);
      }
    }
    list.num_letters = num_letters;
    list.words.push_back(word);
  }

  // Error checking.
  if (list.words.size() == 0 || list.num_letters == 0) {
    die("Empty word list.");
  }

  return list;
}

// Abstract base class for a thing that can play Wordle.
class Strategy {
public:
  // Return the next guess to make.
  virtual std::string MakeGuess() = 0;

  // Process that the given guess got the given response.
  virtual void ProcessResponse(const std::string& guess, const Response& response) = 0;

  virtual ~Strategy() {}
};

class ArbitraryValid : public Strategy {
public:
  ArbitraryValid(const WordList& word_list) : word_list_(word_list) {
    set_ = word_list_.as_set();
  }

  std::string MakeGuess() override {
    if (set_.empty()) {
      die("Can't make a guess if there are no more possible words!");
    }

    // Just arbitrarily pick a word that is still valid.
    const int choice = rand() % set_.size();
    return word_list_.words[set_[choice]];
  }

  void ProcessResponse(const std::string& guess, const Response& response) override {
    // Remove all words that don't conform to the guess.
    set_ = FilterWordSet(word_list_, set_, guess, response);
    // std::cout << "Possible words are now: " << WordSetToString(word_list_, set_) << std::endl;
  }

private:
  // The full list of possible words.
  const WordList& word_list_;

  // The current set of words that are still possible.
  WordSet set_;
};

class MaxEntropy : public Strategy {
public:
  MaxEntropy(const WordList& word_list) : word_list_(word_list) {
    set_ = word_list_.as_set();
  }

  std::string MakeGuess() override {
    if (set_.empty()) {
      die("Can't make a guess if there are no more possible words!");
    }

    // If we know the answer, guess it!
    if (set_.size() == 1) {
      return word_list_.words[set_[0]];
    }

    // Find the max-entropy guess.
    double best_entropy = -1;
    const std::string* best_guess = nullptr;
    for (const std::string& guess : word_list_.words) {
      // Consider each possible remaining word, see what response it would
      // elicit with this guess, and record the distribution over responses.
      std::unordered_map<int, int> response_counts;
      for (const int i : set_) {
	const std::string target = word_list_.words[i];
	const Response response = ScoreGuess(guess, target);
	++response_counts[ResponseToCode(response)];
      }
      // Compute the entropy of this response distribution.
      const double entropy = ComputeEntropy(response_counts);
      // If this entropy is the best so far, record it.
      if (entropy > best_entropy) {
	best_entropy = entropy;
	best_guess = &guess;
      }
    }

    if (best_guess == nullptr) {
      die("All guesses have negative entropy? Bug.");
    }
    return *best_guess;
  }

  void ProcessResponse(const std::string& guess, const Response& response) override {
    // Remove all words that don't conform to the guess.
    set_ = FilterWordSet(word_list_, set_, guess, response);
  }

private:
  // The full list of possible words.
  const WordList& word_list_;

  // The current set of words that are still possible.
  WordSet set_;
};

std::unique_ptr<Strategy> MakeStrategy(const std::string& name,
				       const WordList& word_list) {
  if (name == "ArbitraryValid") {
    return std::make_unique<ArbitraryValid>(word_list);
  } else if (name == "MaxEntropy") {
    return std::make_unique<MaxEntropy>(word_list);
  } else {
    die("Unrecognized strategy name: " + name);
  }
}

int SelfPlay(const std::string& target, Strategy& strategy) {
  std::string guess = "";
  int count = 0;
  while (guess != target) {
    guess = strategy.MakeGuess();
    Response response = ScoreGuess(guess, target);
    DisplayResponse(guess, response);
    strategy.ProcessResponse(guess, response);
    ++count;
  }
  return count;
}

int main(int argc, char* argv[]) {
  std::cout << "Hello world!" << std::endl;
  std::cout << "Wordle is not implemented yet :)" << std::endl;
  const WordList list = ReadWordList("wordlist");
  std::cout << "Number of letters is: " << list.num_letters << std::endl;
  std::cout << "Words are:" << std::endl;
  for (const auto& word : list.words) {
    std::cout << word << std::endl;
  }

  const std::string target = "wince";

  const auto score_and_display = [&](const std::string& guess) {
    const Response response = ScoreGuess(guess, target);
    DisplayResponse(guess, response);
  };

  score_and_display("iotas");
  score_and_display("brief");
  score_and_display("pluck");
  score_and_display("wench");
  score_and_display("wince");

  const WordSet all = list.as_set();
  std::cout << WordSetToString(list, FilterWordSet(list, all, "rxxxx", {PARTIAL_MATCH, NO_MATCH, NO_MATCH, NO_MATCH, NO_MATCH}))
	    << std::endl;

  std::cout << std::endl << std::endl << std::endl;

  for (int i = 0; i < 5; i++) {
    const std::string& word = list.words[rand() % list.words.size()];
    std::cout << "Secret word is: " << word << std::endl;
    std::unique_ptr<Strategy> strategy = MakeStrategy("MaxEntropy", list);
    const int guesses = SelfPlay(word, *strategy);
    std::cout << "Guessed in " << guesses << std::endl;
  }

  return 0;
}
