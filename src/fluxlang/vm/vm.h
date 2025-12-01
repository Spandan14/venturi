#pragma once

namespace Flux {
enum class OpCode {
  PUSH_CONST,
  LOAD_I,
  LOAD_J,
  LOAD_K,
  ADD,
  SUB,
  MUL,
  DIV,
  AND,
  OR,
  LT,
  LE,
  GT,
  GE,
  EQ,
  NEQ,
  SIN,
  COS,
  TAN,
  ABS,
  SQRT,
  LOG,
  EXP,
};

struct Instruction {
  OpCode opcode;
  float operand; // for PUSH_CONST

  Instruction(OpCode op, float oper = 0.0f) : opcode(op), operand(oper) {}
};
} // namespace Flux
