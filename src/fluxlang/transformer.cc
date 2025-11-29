#include <fluxlang/transformer.h>

std::unique_ptr<Flux::Script> FluxASTTransformer::transform() {
  std::vector<std::unique_ptr<Flux::Statement>> statements;

  std::vector<std::shared_ptr<peg::AstBase<peg::Ast::EmptyType>>> node_queue;

  node_queue.push_back(
      std::make_shared<peg::AstBase<peg::Ast::EmptyType>>(peg_root));
  while (!node_queue.empty()) {
    auto current_node = node_queue.back();
    node_queue.pop_back();

    std::cout << "Visiting node: " << current_node->name << std::endl;

    // if (current_node->nodes.empty()) {
    std::cout << "Token: " << current_node->token << std::endl;
    //   continue;
    // }

    for (auto &child : current_node->nodes) {
      node_queue.push_back(child);
    }
  }

  return std::make_unique<Flux::Script>(std::move(statements));
}
