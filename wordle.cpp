#include <fstream>
#include <iostream>
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
void DisplayWordSet(const WordList& word_list, const WordSet& set) {
  for (const int i : set) {
    std::cout << word_list.words[i] << std::endl;
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
  DisplayWordSet(list, FilterWordSet(list, all, "rxxxx", {PARTIAL_MATCH, NO_MATCH, NO_MATCH, NO_MATCH, NO_MATCH}));

  return 0;
}
