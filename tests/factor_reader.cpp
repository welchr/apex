#include "factor_reader.h"
using namespace std;

void FactorReader::load(const string &filepath) {
  if (!filepath_exists(filepath)) {
    throw std::runtime_error("Could not find file: " + filepath);
  }

  unique_ptr<istream> file;
  boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
  ifstream fs(filepath, ios_base::in | ios_base::binary);

  // Check if the file is gzipped by reading the first two bytes
  char bytes[2];
  fs.read(bytes, 2);
  fs.seekg(0, ios_base::beg); // Reset file pointer to the beginning

  if (bytes[0] == '\x1F' && bytes[1] == '\x8B') { // Gzip signature
    inbuf.push(boost::iostreams::gzip_decompressor());
    inbuf.push(fs);
    file = make_unique<istream>(&inbuf);
  } else {
    file = make_unique<ifstream>(move(fs));
  }

  string line;
  auto line_separator = regex("[ \t]");
  auto field_separator = regex(",");

  // Line regexes
  auto regex_header = regex("#id");

  bool header_done = false;
  while (getline(*file, line)) {
    smatch match;
    if (!header_done) {
      if (regex_search(line, match, regex_header)) {
        // Copy all values after the first column (first entry is just "#id")
        auto first_delimiter_pos = line.find_first_of(" \t");
        if (first_delimiter_pos != string::npos) {
          // Start from the character after the first delimiter
          auto start_pos = line.begin() + first_delimiter_pos + 1;
          copy(sregex_token_iterator(start_pos, line.end(), line_separator, -1),
               sregex_token_iterator(),
               back_inserter(samples));
        }
        header_done = true;
      }
    }
    else {
      // Begin parsing rows of covariates
      if (line.substr(0,1) == "#") { continue; }

      vector<string> tokens;
      // Copy tokens to a vector of strings, skipping the first token.
      auto begin_it = sregex_token_iterator(line.begin(), line.end(), line_separator, -1);
      auto end_it = sregex_token_iterator();
      if (begin_it != end_it) {
        // There is at least one token. Skip the first token.
        ++begin_it;
      }
      copy(begin_it, end_it, back_inserter(tokens));

      // Convert the tokens to doubles and store them in another vector.
      vector<double> values;
      transform(tokens.begin(), tokens.end(), back_inserter(values), [](const string& token) {
        return stod(token); // Convert each token to a double.
      });

      // Covariate name
      covariates.push_back(tokens.at(0));

      // Insert
      rows.push_back(values);
    }
  }

  if (rows.empty()) {
    throw std::runtime_error("No results were read from file: " + filepath);
  }
}

FactorReader::FactorReader(const string &file) {
  load(file);
}

bool FactorReader::operator==(const FactorReader& other) {
  if (rows.size() != other.rows.size()) {
    return false;
  }

  // Compare row/column values
  const size_t nrows = rows.size();
  for (size_t i = 0; i < nrows; i++) {
    auto& rec_a = rows[i];
    auto& rec_b = other.rows[i];

    if (rec_a.size() != rec_b.size()) {
      return false;
    }

    for (size_t j = 0; j < rec_a.size(); j++) {
      if (!approx_equal(rec_a[j], rec_b[j])) {
        return false;
      }
    }
  }

  // Compare covariate names
  if (covariates.size() != other.covariates.size()) {
    return false;
  }
  for (size_t i = 0; i < covariates.size(); i++) {
    if (covariates[i] != other.covariates[i]) {
      return false;
    }
  }

  // Compare sample names
  if (samples.size() != other.samples.size()) {
    return false;
  }
  for (size_t i = 0; i < samples.size(); i++) {
    if (samples[i] != other.samples[i]) {
      return false;
    }
  }

  return true;
}