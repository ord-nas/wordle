#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
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
constexpr int MAX_VERBOSE_SET_SIZE = 15;

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

int CorrectGuessCode() {
  Response response;
  response.fill(EXACT_MATCH);
  return ResponseToCode(response);
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
  double entropy = 0.0;
  for (const int entry : distribution) {
    if (entry > 0) {
      const double P = 1.0 * entry / N;
      entropy -= P * log(P);
    }
  }

  return entropy;
}

// The number of expected additional guesses required once we get back a
// response according to this distribution.
double ComputeExpectedGuesses(const ResponseDistribution& distribution) {
  // First count the total number of entries in the distribution.
  int N = 0;
  for (const int entry : distribution) {
    N += entry;
  }

  // Then tally up expected guesses.
  const int correct_guess_code = CorrectGuessCode();
  double guess_sum = 0.0;
  for (int response_code = 0; response_code < distribution.size(); response_code++) {
    const int count = distribution[response_code];
    if (count > 0) {
      const double P = 1.0 * count / N;
      double expected_guesses;
      if (response_code == correct_guess_code) {
	expected_guesses = 0.0;
      } else if (count == 1) {
	expected_guesses = 1.0;
      } else if (count == 2) {
	expected_guesses = 1.5;
      } else {
	// From logarithmic fit on data from `python3 analyze.py stats/max_entropy_counts.csv`
	// expected_guesses = 1.65 + 0.291 * log(count);
	// From logarithmic fit on data from `python3 analyze.py stats/min_expected_guess_counts_0.csv`
	expected_guesses = 1.51 + 0.289 * log(count);
	// Note: `python3 analyze.py stats/min_expected_guess_counts_1.csv`
	// shows basically the same equation, so we have roughly hit a fixed
	// point.
      }
      guess_sum += P * expected_guesses;
    }
  }

  return guess_sum;
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

bool InWordList(const WordList& list, const std::string& word) {
  for (const std::string& entry : list.words) {
    if (entry == word) {
      return true;
    }
  }
  return false;
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

  // Return the number of words that still remain as possibilities given the
  // responses so far.
  int NumRemainingWords() const { return set_.size(); }

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
    const double entropy = ComputeEntropy(distribution);
    const double expected_guesses = ComputeExpectedGuesses(distribution);

    // Now generate the rationale.
    if (set_.size() > MAX_VERBOSE_SET_SIZE) {
      // If the set is too big, just summarize.
      ss << "Words left: " << set_.size() << std::endl;
      ss << "Guess " << guess << " has " << possible_responses << " responses" << std::endl;
      ss << "Entropy = " << entropy << std::endl;
      ss << "Expected remaining guesses = " << expected_guesses;
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
      ss << "}" << std::endl;
      ss << "Entropy = " << entropy << std::endl;
      ss << "Expected remaining guesses = " << expected_guesses;
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

class MinExpectedGuesses : public BestResponseDistribution {
public:
  MinExpectedGuesses(const WordList& word_list) : BestResponseDistribution(word_list) {}

  bool IsMaximizer() const override { return false; }
  double ScoreDistribution(const ResponseDistribution& distribution) const override {
    return ComputeExpectedGuesses(distribution);
  }
};

std::unique_ptr<Strategy> MakeStrategy(const std::string& name,
				       const WordList& word_list) {
  if (name == "ArbitraryValid") {
    return std::make_unique<ArbitraryValid>(word_list);
  } else if (name == "MaxEntropy") {
    return std::make_unique<MaxEntropy>(word_list);
  } else if (name == "MinExpectedGuesses") {
    return std::make_unique<MinExpectedGuesses>(word_list);
  } else {
    die("Unrecognized strategy name: " + name);
  }
}

enum Verbosity {
  SILENT = 0,
  NORMAL = 1,
  VERBOSE = 2,
};

Verbosity ToVerbosity(const std::string& str) {
  if (str == "SILENT") {
    return SILENT;
  } else if (str == "NORMAL") {
    return NORMAL;
  } else if (str == "VERBOSE") {
    return VERBOSE;
  } else {
    die("Unrecognized verbosity: " + str);
  }
}

struct GameOutcome {
  int guess_count = 0;
  // remaining_word_history[i] means the number of remaining words after i
  // receiving i responses.
  std::vector<int> remaining_word_history;
};

GameOutcome SelfPlay(const std::string& target,
	     Strategy& strategy,
	     std::vector<std::string> forced_guesses,
	     Verbosity verbosity) {
  Guess guess;
  GameOutcome outcome;
  while (guess.word != target) {
    outcome.remaining_word_history.push_back(strategy.NumRemainingWords());
    if (!forced_guesses.empty()) {
      guess.word = forced_guesses[0];
      guess.reasoning = "Forced guess";
      forced_guesses.erase(forced_guesses.begin());
    } else {
      guess = strategy.MakeGuess();
    }
    Response response = ScoreGuess(guess.word, target);
    if (verbosity >= VERBOSE && !guess.reasoning.empty()) {
      std::cout << guess.reasoning << std::endl;
    }
    if (verbosity >= NORMAL) {
      std::cout << ColorGuess(guess.word, response) << std::endl;
    }
    strategy.ProcessResponse(guess.word, response);
    ++outcome.guess_count;
  }
  return outcome;
}

void SelfPlayLoop(const WordList& list, const Flags& flags) {
  const std::string strategy_name = flags.Get("strategy", /*default=*/"MaxEntropy");
  const std::vector<std::string> forced_guesses = Split(flags.Get("forced_guesses", /*default=*/""), ',');
  const Verbosity verbosity = ToVerbosity(flags.Get("verbosity", "NORMAL"));

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
    const GameOutcome outcome = SelfPlay(word, *strategy, forced_guesses, verbosity);
    std::cout << "Guessed in " << outcome.guess_count << std::endl;
  }
}

int HumanPlay(const WordList& list, const std::string& target) {
  int count = 0;
  std::string guess;
  while (guess != target) {
    // Get word.
    std::cout << "Enter guess: ";
    std::cin >> guess;
    if (!InWordList(list, guess)) {
      std::cout << "Not in word list!" << std::endl;
      continue;
    }

    // Show response.
    Response response = ScoreGuess(guess, target);
    std::cout << ColorGuess(guess, response) << std::endl;
    ++count;
  }
  return count;
}

void HumanPlayLoop(const WordList& list) {
  while (true) {
    std::string word = list.words[rand() % list.words.size()];
    const int guesses = HumanPlay(list, word);
    std::cout << "Guessed in " << guesses << std::endl;
  }
}

std::vector<std::string> Choose(const WordList& list, int N) {
  if (N >= list.words.size()) {
    return list.words;
  }

  std::vector<std::string> remaining = list.words;
  std::vector<std::string> selected;
  for (int i = 0; i < N; i++) {
    const int choice = rand() % remaining.size();
    selected.push_back(remaining[choice]);
    remaining[choice] = remaining.back();
    remaining.pop_back();
  }

  return selected;
}

class ProgressReporter {
public:
  ProgressReporter(Verbosity verbosity, int num_rounds)
    : verbosity_(verbosity), num_rounds_(num_rounds) {}

  void Report(int round, const std::string& word, const std::string& strategy) {
    const std::string message = GetMessage(round, word, strategy);
    if (verbosity_ >= NORMAL) {
      std::cout << message << std::endl;
    } else {
      for (int i = 0; i < last_message_.size(); i++) {
	std::cout << "\b";
      }
      std::cout << message << std::flush;
      last_message_ = message;
    }
  }

  void Done() {
    if (verbosity_ < NORMAL) {
      std::cout << std::endl;
    }
  }

private:
  std::string GetMessage(int round, const std::string& word, const std::string& strategy) {
    std::ostringstream ss;
    ss << "[" << (round + 1) << "/" << num_rounds_ << "]: " << word << " - " << strategy;
    return ss.str();
  }

  Verbosity verbosity_;
  int num_rounds_;
  std::string last_message_;
};

struct StrategyStats {
  std::vector<int> guess_count_history;
  std::map<int, std::vector<int>> remaining_words_to_guesses;
};

void FinalizeStats(std::vector<StrategyStats>& stats_list) {
  // Sort all the guess lists in remaining_words_to_guesses.
  for (auto& stats : stats_list) {
    for (auto& entry : stats.remaining_words_to_guesses) {
      std::sort(entry.second.begin(), entry.second.end());
    }
  }
}

void CollectStats(const WordList& list, const Flags& flags) {
  // Pull in some flag values.
  const std::vector<std::string> strategy_names = Split(flags.Get("strategies"), ',');
  const std::vector<std::string> forced_guesses = Split(flags.Get("forced_guesses", /*default=*/""), ',');
  const Verbosity verbosity = ToVerbosity(flags.Get("verbosity", "SILENT"));
  const std::string filename = flags.Get("out_file", "");
  int rounds = flags.GetInt("rounds", "100");

  // Choose some words to play.
  if (rounds <= 0) {
    rounds = list.words.size();
  }
  std::vector<std::string> words = Choose(list, rounds);
  rounds = words.size();

  // Play each of the words through each of the strategies.
  ProgressReporter progress(verbosity, rounds);
  std::vector<StrategyStats> stats(strategy_names.size());
  for (int i = 0; i < rounds; i++) {
    const std::string& word = words[i];
    for (int j = 0; j < strategy_names.size(); j++) {
      progress.Report(i, word, strategy_names[j]);
      std::unique_ptr<Strategy> strategy = MakeStrategy(strategy_names[j], list);
      const GameOutcome outcome = SelfPlay(word, *strategy, forced_guesses, verbosity);
      // Add the outcome to overall stats.
      stats[j].guess_count_history.push_back(outcome.guess_count);
      for (int i = 0; i < outcome.remaining_word_history.size(); i++) {
	const int remaining_words = outcome.remaining_word_history[i];
	const int remaining_guesses = outcome.guess_count - i;
	stats[j].remaining_words_to_guesses[remaining_words].push_back(remaining_guesses);
      }
    }
  }
  FinalizeStats(stats);
  progress.Done();
  std::cout << std::endl << std::endl;

  // Compile a report for stdout.
  for (int j = 0; j < strategy_names.size(); j++) {
    std::cout << strategy_names[j] << std::endl;
    int total_guesses = 0;
    std::map<int, int> guess_distribution;
    for (const int guess_count : stats[j].guess_count_history) {
      ++guess_distribution[guess_count];
      total_guesses += guess_count;
    }
    std::cout << "Average guesses: " << (1.0 * total_guesses / rounds) << std::endl;
    for (const auto& entry : guess_distribution) {
      std::cout << entry.first << " guesses: " << (100.0 * entry.second / rounds) << " %" << std::endl;
    }
    std::cout << std::endl;
  }

  // Compile a report for file.
  if (!filename.empty()) {
    std::ofstream file(filename);
    file << "strategy,remaining_words,guess_count" << std::endl;
    for (int j = 0; j < strategy_names.size(); j++) {
      for (const auto& entry : stats[j].remaining_words_to_guesses) {
	for (const int guess_count : entry.second) {
	  file << strategy_names[j] << "," << entry.first << "," << guess_count << std::endl;
	}
      }
    }
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
  } else if (mode == "human_play") {
    HumanPlayLoop(list);
  } else if (mode == "stats") {
    CollectStats(list, flags);
  } else if (mode == "ai_play") {
    die("Unimplemented.");
  } else {
    die("Unrecognized mode: " + mode);
  }

  return 0;
}
