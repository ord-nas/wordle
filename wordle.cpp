#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

struct WordList {
  std::vector<std::string> words;
  int num_letters = 0;
};

enum Outcome {
  EXACT_MATCH,
  PARTIAL_MATCH,
  NO_MATCH,
};

using Response = std::vector<Outcome>;

void die(const std::string& message) {
  std::cout << message << std::endl;
  exit(1);
}

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

  return 0;
}
