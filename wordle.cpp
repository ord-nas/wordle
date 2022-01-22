#include <fstream>
#include <iostream>
#include <limits>
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

constexpr int NUM_LETTERS = 5;
constexpr int NUM_RESPONSES = 243; // 3 ^ NUM_LETTERS

// Guess outcomes for a full word. Response[i] is the outcome for guess[i].
using Response = std::array<Outcome, NUM_LETTERS>;

// Print the given error message and abort.
void die(const std::string& message) {
  std::cout << message << std::endl;
  exit(1);
}

// Score guess against target.
Response ScoreGuess(const std::string& guess, std::string target) {
  Response response;
  for (int i = 0; i < NUM_LETTERS; i++) {
    if (guess[i] == target[i]) {
      response[i] = EXACT_MATCH;
    } else {
      response[i] = NO_MATCH;
      // Try to find a partial match.
      for (int j = i + 1; j < NUM_LETTERS; j++) {
	if (guess[i] == target[j] && guess[j] != target[j]) {
	  target[j] = '.';
	  response[i] = PARTIAL_MATCH;
	  break;
	}
      }
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

using ResponseDistribution = std::array<int, NUM_RESPONSES>;

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

double ComputeEntropy(const ResponseDistribution& distribution) {
  // First count the total number of entries in the distribution.
  int N = 0;
  for (const int entry : distribution) {
    N += entry;
  }

  // Then tally up entropy.
  double entropy;
  for (const int entry : distribution) {
    if (entry > 0) {
      const double P = 1.0 * entry / N;
      entropy -= P * log(P);
    }
  }

  return entropy;
}

std::string ColorGuess(const std::string& guess, const Response& response) {
  if (guess.size() != response.size()) {
    die("Can't display guess and target of different sizes. guess=" + guess);
  }

  std::ostringstream ss;
  for (int i = 0; i < guess.size(); i++) {
    switch (response[i]) {
      case NO_MATCH:
	ss << "\033[37;40m" << guess[i] << "\033[0m";
	break;
      case EXACT_MATCH:
	ss << "\033[37;42m" << guess[i] << "\033[0m";
	break;
      case PARTIAL_MATCH:
	ss << "\033[37;43m" << guess[i] << "\033[0m";
	break;
    }
  }

  return ss.str();
}

// Display the response for this guess, using ANSI color codes.
void DisplayResponse(const std::string& guess, const Response& response) {
  std::cout << ColorGuess(guess, response) << std::endl;
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
    if (word.size() != NUM_LETTERS) {
      die("Got word with wrong number of letters: " + word);
    }
    for (char c : word) {
      if (c < 'a' || c > 'z') {
	die("Got invalid character in word: " + word);
      }
    }
    list.words.push_back(word);
  }

  // Error checking.
  if (list.words.size() == 0) {
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
      ResponseDistribution distribution;
      distribution.fill(0);
      for (const int i : set_) {
	const std::string& target = word_list_.words[i];
	const Response response = ScoreGuess(guess, target);
	++distribution[ResponseToCode(response)];
      }
      // Compute the entropy of this response distribution.
      const double entropy = ComputeEntropy(distribution);
      // If this entropy is the best so far, record it.
      if (entropy > best_entropy) {
      	best_entropy = entropy;
      	best_guess = &guess;
      }
    }

    if (best_guess == nullptr) {
      die("All guesses have negative entropy? Bug.");
    }

    if (set_.size() <= 10) {
      std::unordered_map<std::string, std::vector<std::string>> response_to_targets;
      for (const int i : set_) {
	const std::string& target = word_list_.words[i];
	const Response response = ScoreGuess(*best_guess, target);
	response_to_targets[ColorGuess(*best_guess, response)].push_back(target);
      }

      std::cout << "For guess " << *best_guess << ", potential reponses are: {";
      bool outer_first = true;
      for (const auto& entry : response_to_targets) {
	if (!outer_first) std::cout << ", ";
	std::cout << entry.first << " => (";
	bool inner_first = true;
	for (const auto& target : entry.second) {
	  if (!inner_first) std::cout << ", ";
	  std::cout << target;
	  inner_first = false;
	}
	std::cout << ")";
	outer_first = false;
      }
      std::cout << "}" << std::endl;
    }

    return *best_guess;
  }

  void ProcessResponse(const std::string& guess, const Response& response) override {
    // Remove all words that don't conform to the guess.
    set_ = FilterWordSet(word_list_, set_, guess, response);
    if (set_.size() > 10) {
      std::cout << "Words left: " << set_.size() << std::endl;
    } else {
      std::cout << "Words left: " << set_.size() << " " << WordSetToString(word_list_, set_) << std::endl;
    }
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
  const WordList list = ReadWordList("wordlist");
  for (const std::string word : {"wince", "prick", "robot", "point", "proxy", "shire", "solar", "panic", "tangy"}) {
    // const std::string& word = list.words[rand() % list.words.size()];
    std::cout << "Secret word is: " << word << std::endl;
    std::unique_ptr<Strategy> strategy = MakeStrategy("MaxEntropy", list);
    const int guesses = SelfPlay(word, *strategy);
    std::cout << "Guessed in " << guesses << std::endl;
  }
  return 0;
}
