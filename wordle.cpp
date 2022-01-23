#include <fstream>
#include <iostream>
#include <limits>
#include <math.h>
#include <memory>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <unordered_map>
#include <vector>

// Print the given error message and abort.
void die(const std::string& message) {
  std::cout << message << std::endl;
  exit(1);
}

// Split a string on the given delimiter character.
std::vector<std::string> Split(const std::string& str, char delimiter) {
  if (str.empty()) {
    return {};
  }

  std::vector<std::string> elements;
  std::string current;
  for (const char c : str) {
    if (c == delimiter) {
      elements.push_back(current);
      current = "";
    } else {
      current += c;
    }
  }
  elements.push_back(current);
  return elements;
}

int ToInt(const std::string& str) {
  int x = 0;
  std::size_t chars_processed = 0;
  try {
    x = std::stoi(str, &chars_processed);
  } catch (...) {
    die("Couldn't convert to int: " + str);
  }

  if (chars_processed != str.size()) {
    die("Couldn't convert to int: " + str);
  }

  return x;
}

class Flags {
public:
  Flags(int argc, char* argv[]) {
    for (int i = 1; i < argc; i++) {
      std::string flag(argv[i]);

      // Flag should start with --
      if (flag.size() <= 2 || flag[0] != '-' || flag[1] != '-') {
	die("Malformed flag: " + flag);
      }

      // Remove the --
      flag.erase(0, 2);

      // Flag should now have form key=value
      std::vector<std::string> parts = Split(flag, '=');
      if (parts.size() != 2 || parts[0].empty() || parts[1].empty()) {
	die("Malformed flag parts: " + flag);
      }

      // Store results.
      values_[parts[0]] = parts[1];
    }
  }

  // Check if flags was specified.
  bool Has(const std::string& key) const {
    return values_.count(key) > 0;
  }

  // String accessors.
  std::string Get(const std::string& key) const {
    const auto it = values_.find(key);
    if (it == values_.end()) {
      die("Missing required flag: " + key);
    }
    return it->second;
  }
  std::string Get(const std::string& key, const std::string& default_value) const {
    const auto it = values_.find(key);
    if (it == values_.end()) {
      return default_value;
    }
    return it->second;
  }

  // Int accessors.
  int GetInt(const std::string& key) const {
    return ToInt(Get(key));
  }
  int GetInt(const std::string& key, const std::string& default_value) const {
    return ToInt(Get(key, default_value));
  }

  // Print all the flag values to stdout.
  void Print() const {
    std::cout << "Flags {" << std::endl;
    for (const auto& entry : values_) {
      std::cout << "    " << entry.first << " => " << entry.second << std::endl;
    }
    std::cout << "}" << std::endl;
  }

private:
  std::unordered_map<std::string, std::string> values_;
};

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

// Score guess against target.
Response ScoreGuess(const std::string& guess, std::string target) {
  Response response;
  for (int i = 0; i < NUM_LETTERS; i++) {
    if (guess[i] == target[i]) {
      response[i] = EXACT_MATCH;
    } else {
      response[i] = NO_MATCH;
      // Try to find a partial match.
      for (int j = 0; j < NUM_LETTERS; j++) {
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

double CountEntries(const ResponseDistribution& distribution) {
  // Count how many distinct responses.
  int N = 0;
  for (const int entry : distribution) {
    if (entry > 0) ++N;
  }
  return N;
}

// Color the guess based on the response, using ANSI color codes.
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

void ValidateWord(const std::string& word) {
  if (word.size() != NUM_LETTERS) {
    die("Got word with wrong number of letters: " + word);
  }
  for (char c : word) {
    if (c < 'a' || c > 'z') {
      die("Got invalid character in word: " + word);
    }
  }
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
    ValidateWord(word);
    list.words.push_back(word);
  }

  // Error checking.
  if (list.words.size() == 0) {
    die("Empty word list.");
  }

  return list;
}

struct Guess {
  std::string word;
  std::string reasoning;
};

// Base class for a thing that can play Wordle.
class Strategy {
public:
  Strategy(const WordList& word_list) : word_list_(word_list) {
    set_ = word_list_.as_set();
  }

  // Return the next guess to make.
  virtual Guess MakeGuess() = 0;

  // Process that the given guess got the given response.
  virtual void ProcessResponse(const std::string& guess, const Response& response) {
    // Remove all words that don't conform to the guess.
    set_ = FilterWordSet(word_list_, set_, guess, response);
  }

  virtual ~Strategy() {}

protected:
  // The full list of possible words.
  const WordList& word_list_;

  // The current set of words that are still possible.
  WordSet set_;

};

class ArbitraryValid : public Strategy {
public:
  ArbitraryValid(const WordList& word_list) : Strategy(word_list) {}

  Guess MakeGuess() override {
    if (set_.empty()) {
      die("Can't make a guess if there are no more possible words!");
    }

    // Just arbitrarily pick a word that is still valid.
    Guess guess;
    const int choice = rand() % set_.size();
    guess.word =  word_list_.words[set_[choice]];
    return guess;
  }
};

class BestResponseDistribution : public Strategy {
public:
  BestResponseDistribution(const WordList& word_list) : Strategy(word_list) {}

  virtual bool IsMaximizer() const = 0;
  virtual double ScoreDistribution(const ResponseDistribution& distribution) const = 0;

  Guess MakeGuess() override {
    if (set_.empty()) {
      die("Can't make a guess if there are no more possible words!");
    }

    // If we know the answer, guess it!
    if (set_.size() == 1) {
      Guess guess;
      guess.word = word_list_.words[set_[0]];
      guess.reasoning = "Only one word remaining";
      return guess;
    }

    // // Hack for the first round.
    // if (set_.size() == word_list_.words.size()) {
    //   Guess guess;
    //   guess.word = "rates";
    //   guess.reasoning = "Forced guess";
    //   return guess;
    // }

    // Find the best guess.
    double best_score = IsMaximizer() ?
      -std::numeric_limits<double>::infinity() :
      std::numeric_limits<double>::infinity();
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
      // Compute some score over the distribution.
      const double score = ScoreDistribution(distribution);
      // If this score is the best so far, record it.
      if (IsMaximizer() ? (score > best_score) : (score < best_score)) {
      	best_score = score;
      	best_guess = &guess;
      }
    }

    if (best_guess == nullptr) {
      die("All guesses have invalid scores? Bug.");
    }

    Guess guess;
    guess.word = *best_guess;
    guess.reasoning = ExplainGuess(*best_guess);
    return guess;
  }

private:
  std::string ExplainGuess(const std::string& guess) {
    std::ostringstream ss;

    // Count up the number of distinct responses that are possible given this
    // guess.
    ResponseDistribution distribution;
    distribution.fill(0);
    for (const int i : set_) {
      const std::string& target = word_list_.words[i];
      const Response response = ScoreGuess(guess, target);
      ++distribution[ResponseToCode(response)];
    }
    const int possible_responses = CountEntries(distribution);

    // Now generate the rationale.
    if (set_.size() > 10) {
      // If the set is too big, just summarize.
      ss << "Words left: " << set_.size() << std::endl;
      ss << "Guess " << guess << " has " << possible_responses << " responses";
    } else {
      // Otherwise, go into more detail.
      ss << "Words left: " << set_.size() << " " << WordSetToString(word_list_, set_) << std::endl;

      // Partition the set of possible words based on the response to our guess.
      std::unordered_map<std::string, std::vector<std::string>> response_to_targets;
      for (const int i : set_) {
	const std::string& target = word_list_.words[i];
	const Response response = ScoreGuess(guess, target);
	response_to_targets[ColorGuess(guess, response)].push_back(target);
      }

      // Format.
      ss << "Guess " << guess << " has " << possible_responses << " responses {";
      bool outer_first = true;
      for (const auto& entry : response_to_targets) {
	if (!outer_first) ss << ", ";
	ss << entry.first << " => (";
	bool inner_first = true;
	for (const auto& target : entry.second) {
	  if (!inner_first) ss << ", ";
	  ss << target;
	  inner_first = false;
	}
	ss << ")";
	outer_first = false;
      }
      ss << "}";
    }

    return ss.str();
  }
};

class MaxEntropy : public BestResponseDistribution {
public:
  MaxEntropy(const WordList& word_list) : BestResponseDistribution(word_list) {}

  bool IsMaximizer() const override { return true; }
  double ScoreDistribution(const ResponseDistribution& distribution) const override {
    return ComputeEntropy(distribution);
  }
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

enum DisplayMode {
  SILENT = 0,
  NORMAL = 1,
  VERBOSE = 2,
};

int SelfPlay(const std::string& target,
	     Strategy& strategy,
	     std::vector<std::string> forced_guesses,
	     DisplayMode display_mode) {
  Guess guess;
  guess.word = "";
  guess.reasoning = "";
  int count = 0;
  while (guess.word != target) {
    if (!forced_guesses.empty()) {
      guess.word = forced_guesses[0];
      guess.reasoning = "Forced guess";
      forced_guesses.erase(forced_guesses.begin());
    } else {
      guess = strategy.MakeGuess();
    }
    Response response = ScoreGuess(guess.word, target);
    if (display_mode >= VERBOSE && !guess.reasoning.empty()) {
      std::cout << guess.reasoning << std::endl;
    }
    if (display_mode >= NORMAL) {
      std::cout << ColorGuess(guess.word, response) << std::endl;
    }
    strategy.ProcessResponse(guess.word, response);
    ++count;
  }
  return count;
}

void SelfPlayLoop(const WordList& list, const Flags& flags) {
  const std::string strategy_name = flags.Get("strategy", /*default=*/"MaxEntropy");

  const std::vector<std::string> forced_guesses = Split(flags.Get("forced_guesses", /*default=*/""), ',');

  while (true) {
    std::cout << "Enter a word to play, or empty string to pick a random word: ";
    std::string word;
    std::getline(std::cin, word);
    if (word.empty()) {
      word = list.words[rand() % list.words.size()];
      std::cout << "Secret word is: " << word << std::endl;
    } else {
      ValidateWord(word);
    }
    std::unique_ptr<Strategy> strategy = MakeStrategy(strategy_name, list);
    const int guesses = SelfPlay(word, *strategy, forced_guesses, /*display_mode=*/VERBOSE);
    std::cout << "Guessed in " << guesses << std::endl;
  }
}

int main(int argc, char* argv[]) {
  const Flags flags(argc, argv);
  flags.Print();

  // Read the word list.
  const WordList list = ReadWordList(flags.Get("wordlist", /*default=*/"wordlist"));

  // Initialize RNG.
  if (flags.Has("seed")) {
    srand(flags.GetInt("seed"));
  } else {
    srand(time(nullptr));
  }

  // Get the mode.
  const std::string mode = flags.Get("mode", /*default=*/"self_play");
  if (mode == "self_play") {
    SelfPlayLoop(list, flags);
  }

  return 0;
}
