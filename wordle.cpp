#include <fstream>
#include <iostream>
#include <string>
#include <vector>

struct WordList {
  std::vector<std::string> words;
  int num_letters = 0;
};

void die(const std::string& message) {
  std::cout << message << std::endl;
  exit(1);
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
  return 0;
}
