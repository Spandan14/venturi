#include "parser.h"

peg::parser Parser::parser = Parser::_load_grammar(FLUX_LANG_GRAMMAR_PATH);

std::shared_ptr<peg::Ast> Parser::parse(const char *filename) {
  source = load_text_file(filename);

  if (!parser.parse(source.c_str(), ast)) {
    throw std::runtime_error("Failed to parse source file!");
  }

  return ast;
}

peg::parser Parser::_load_grammar(const std::string &grammar_path) {
  std::string grammar = load_text_file(grammar_path);

  peg::parser p(grammar.c_str());

  if (!p) {
    throw std::runtime_error("Failed to build parser from grammar file!");
  }

  p.set_logger([](size_t line, size_t col, const std::string &msg,
                  const std::string &rule) {
    std::cerr << line << ":" << col << ": " << msg << "\n";
  });
  p.enable_ast();

  return p;
}
