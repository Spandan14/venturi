#include "engine/sim3d.h"
#include "fluxlang/ast/nodes.h"
#include "renderer/3d/renderer_3d.h"
#include <fluxlang/runtime.h>
#include <iostream>

void Runtime::run() {
  eval(*_script);

  _prepare_gen();

  if (_config.dim == Flux::DimType::Two) {
    _prepare_2d();
  } else if (_config.dim == Flux::DimType::Three) {
    _prepare_3d();
  } else {
    throw std::runtime_error(
        "Dimension type must be specified before running the simulation!");
  }

  auto last_frame = std::chrono::high_resolution_clock::now();
  while (_config.renderer->should_draw()) {
    auto this_frame = std::chrono::high_resolution_clock::now();
    auto step_ms =
        std::chrono::duration<float>(this_frame - last_frame).count();

    _step(step_ms);
    last_frame = this_frame;

    _config.renderer->render();
    _config.renderer->poll();
  }
}

void Runtime::_prepare_gen() {
  // move 'all' to be at the start of the set, so that it is evaluated first
  // auto it =
  //     std::find(_config.set_order.begin(), _config.set_order.end(), "all");
  // if (it != _config.set_order.end()) {
  //   _config.set_order.erase(it);
  _config.set_order.insert(_config.set_order.begin(), "all");
  // }

  // make membership for all as well
  if (_config.dim == Flux::DimType::Two) {
    auto all_2d = [](int i, int j, float t) { return true; };
    Membership2D all_membership = Membership2D(all_2d);
    _config.set_memberships["all"] = Membership2D(all_membership);
  } else if (_config.dim == Flux::DimType::Three) {
    auto all_3d = [](int i, int j, int k, float t) { return true; };
    Membership3D all_membership = Membership3D(all_3d);
    _config.set_memberships["all"] = Membership3D(all_membership);
  }
}

void Runtime::_prepare_2d() {
  auto sets_lambda = [this](int i, int j, float t) {
    std::vector<std::string> sets;
    for (const auto &set_name : _config.set_order) {
      const auto &membership = _config.set_memberships[set_name];
      if (std::holds_alternative<Membership2D>(membership)) {
        auto &mem_func = std::get<Membership2D>(membership);
        if (mem_func(i, j, t)) {
          sets.push_back(set_name);
        }
      }
    }
    return sets;
  };

  auto density_initializer = [this, sets_lambda](int i, int j, float t) {
    float total_density = 0.0f;

    std::vector<std::string> sets_to_check = sets_lambda(i, j, t);

    for (const auto &set_name : sets_to_check) {
      if (_config.density_values.count(set_name) > 0) {
        float density =
            extract_float_2d(_config.density_values[set_name],
                             static_cast<float>(i), static_cast<float>(j), t);

        total_density = density;
      }
    }

    return total_density;
  };

  auto force_initializer = [this, sets_lambda](int i, int j, float t) {
    vec2d total_force(0.0f, 0.0f);

    std::vector<std::string> sets_to_check = sets_lambda(i, j, t);
    for (const auto &set_name : sets_to_check) {
      if (_config.force_values.count(set_name) > 0) {
        auto &force_values = _config.force_values[set_name]; // vector of Value

        if (force_values->size() != 2) {
          throw std::runtime_error(
              "Force value vector must have exactly 2 components for 2D "
              "simulations!");
        }

        float fx = extract_float_2d((*force_values)[0], static_cast<float>(i),
                                    static_cast<float>(j), t);
        float fy = extract_float_2d((*force_values)[1], static_cast<float>(i),
                                    static_cast<float>(j), t);

        total_force = vec2d(fx, fy);
      }
    }

    return total_force;
  };

  auto solid_initializer = [this, sets_lambda](int i, int j, float t) {
    bool is_solid = false;

    std::vector<std::string> sets_to_check = sets_lambda(i, j, t);
    for (const auto &set_name : sets_to_check) {
      if (_config.solid_values.count(set_name) > 0) {
        auto &solid_value = _config.solid_values[set_name];
        float val = extract_float_2d(solid_value, static_cast<float>(i),
                                     static_cast<float>(j), t);

        if (val != 0.0f) {
          is_solid = true;
          break;
        }
      }
    }

    return is_solid;
  };

  auto flow_generator = [this, sets_lambda](int i, int j, float t) {
    float total_flow = 0.0f;

    std::vector<std::string> sets_to_check = sets_lambda(i, j, t);
    for (const auto &set_name : sets_to_check) {
      if (_config.flow_values.count(set_name) > 0) {
        auto &flow_value = _config.flow_values[set_name];
        float fx = extract_float_2d(flow_value, static_cast<float>(i),
                                    static_cast<float>(j), t);

        total_flow = fx;
      }
    }

    return total_flow;
  };

  auto flow_ratio_generator = [this, sets_lambda](int i, int j, float t) {
    float total_flow_ratio = 1.0f;

    std::vector<std::string> sets_to_check = sets_lambda(i, j, t);
    for (const auto &set_name : sets_to_check) {
      if (_config.flow_ratio_values.count(set_name) > 0) {
        auto &flow_ratio_value = _config.flow_ratio_values[set_name];
        float fx = extract_float_2d(flow_ratio_value, static_cast<float>(i),
                                    static_cast<float>(j), t);

        total_flow_ratio = fx;
      }
    }

    return total_flow_ratio;
  };

  auto &sim = std::get<std::unique_ptr<Simulation<2>>>(_config.simulation);
  sim->initialize_density(density_initializer);
  sim->initialize_forces(force_initializer);
  sim->initialize_solids(solid_initializer);
  sim->initialize_flows(flow_generator);
  sim->initialize_flow_ratios(flow_ratio_generator);
}

void Runtime::_prepare_3d() {

  auto sets_lambda = [this](int i, int j, int k, float t) {
    std::vector<std::string> sets;
    for (const auto &set_name : _config.set_order) {
      const auto &membership = _config.set_memberships[set_name];
      if (std::holds_alternative<Membership3D>(membership)) {
        auto &mem_func = std::get<Membership3D>(membership);
        if (mem_func(i, j, k, t)) {
          sets.push_back(set_name);
        }
      }
    }
    return sets;
  };

  auto density_initializer = [this, sets_lambda](int i, int j, int k, float t) {
    float total_density = 0.0f;

    std::vector<std::string> sets_to_check = sets_lambda(i, j, k, t);
    for (const auto &set_name : sets_to_check) {
      if (_config.density_values.count(set_name) > 0) {
        float density = extract_float_3d(
            _config.density_values[set_name], static_cast<float>(i),
            static_cast<float>(j), static_cast<float>(k), t);

        total_density = density;
      }
    }

    return total_density;
  };

  auto force_initializer = [this, sets_lambda](int i, int j, int k, float t) {
    vec3d total_force(0.0f, 0.0f, 0.0f);

    std::vector<std::string> sets_to_check = sets_lambda(i, j, k, t);
    for (const auto &set_name : sets_to_check) {
      if (_config.force_values.count(set_name) > 0) {
        auto &force_values = _config.force_values[set_name]; // vector of Value

        if (force_values->size() != 3) {
          throw std::runtime_error(
              "Force value vector must have exactly 3 components for 3D "
              "simulations!");
        }

        float fx =
            extract_float_3d((*force_values)[0], static_cast<float>(i),
                             static_cast<float>(j), static_cast<float>(k), t);
        float fy =
            extract_float_3d((*force_values)[1], static_cast<float>(i),
                             static_cast<float>(j), static_cast<float>(k), t);
        float fz =
            extract_float_3d((*force_values)[2], static_cast<float>(i),
                             static_cast<float>(j), static_cast<float>(k), t);

        total_force = vec3d(fx, fy, fz);
      }
    }

    return total_force;
  };

  auto solid_initializer = [this, sets_lambda](int i, int j, int k, float t) {
    bool is_solid = false;

    std::vector<std::string> sets_to_check = sets_lambda(i, j, k, t);
    for (const auto &set_name : sets_to_check) {
      if (_config.solid_values.count(set_name) > 0) {
        auto &solid_value = _config.solid_values[set_name];
        float val =
            extract_float_3d(solid_value, static_cast<float>(i),
                             static_cast<float>(j), static_cast<float>(k), t);

        if (val != 0.0f) {
          is_solid = true;
          break;
        }
      }
    }

    return is_solid;
  };

  auto flow_generator = [this, sets_lambda](int i, int j, int k, float t) {
    float total_flow = 0.0f;

    std::vector<std::string> sets_to_check = sets_lambda(i, j, k, t);
    for (const auto &set_name : sets_to_check) {
      if (_config.flow_values.count(set_name) > 0) {
        auto &flow_value = _config.flow_values[set_name];
        float fx =
            extract_float_3d(flow_value, static_cast<float>(i),
                             static_cast<float>(j), static_cast<float>(k), t);

        total_flow = fx;
      }
    }

    return total_flow;
  };

  auto flow_ratio_generator = [this, sets_lambda](int i, int j, int k,
                                                  float t) {
    float total_flow_ratio = 1.0f;

    std::vector<std::string> sets_to_check = sets_lambda(i, j, k, t);
    for (const auto &set_name : sets_to_check) {
      if (_config.flow_ratio_values.count(set_name) > 0) {
        auto &flow_ratio_value = _config.flow_ratio_values[set_name];
        float fx =
            extract_float_3d(flow_ratio_value, static_cast<float>(i),
                             static_cast<float>(j), static_cast<float>(k), t);

        total_flow_ratio = fx;
      }
    }

    return total_flow_ratio;
  };

  auto &sim = std::get<std::unique_ptr<Simulation<3>>>(_config.simulation);
  sim->initialize_density(density_initializer);
  sim->initialize_forces(force_initializer);
  sim->initialize_solids(solid_initializer);
  sim->initialize_flows(flow_generator);
  sim->initialize_flow_ratios(flow_ratio_generator);
}

void Runtime::_step(float dt) {
  if (std::holds_alternative<std::unique_ptr<Simulation<2>>>(
          _config.simulation)) {
    auto &sim = std::get<std::unique_ptr<Simulation<2>>>(_config.simulation);
    sim->step(dt);
  } else if (std::holds_alternative<std::unique_ptr<Simulation<3>>>(
                 _config.simulation)) {
    auto &sim = std::get<std::unique_ptr<Simulation<3>>>(_config.simulation);
    sim->step(dt);
  }
}

void Runtime::eval(Flux::Script &node) {
  for (auto &stmt : node.statements) {
    eval(*stmt);
  }
};

std::function<float(const UnifiedValue &, float, float, float)>
    Runtime::extract_float_2d =
        [](const UnifiedValue &binding, float i, float j, float t) -> float {
  if (std::holds_alternative<Value>(binding)) {
    const auto &val = std::get<Value>(binding);
    return std::get<Value2D>(val)(i, j, t);
  }
  throw std::runtime_error("Invalid binding type in RangeExpression!");
};

std::function<float(const UnifiedValue &, float, float, float, float)>
    Runtime::extract_float_3d = [](const UnifiedValue &binding, float i,
                                   float j, float k, float t) -> float {
  if (std::holds_alternative<Value>(binding)) {
    // Determine if it's 2D or 3D based on what's stored
    const auto &val = std::get<Value>(binding);
    return std::get<Value3D>(val)(i, j, k, t);
  }
  throw std::runtime_error("Invalid binding type in RangeExpression!");
};

UnifiedValue Runtime::eval(Flux::Expression &node) {
  if (auto literal = dynamic_cast<Flux::LiteralExpression *>(&node)) {
    return eval(*literal);
  }

  if (auto gen_var = dynamic_cast<Flux::GenVariableExpression *>(&node)) {
    return eval(*gen_var);
  }

  if (auto range = dynamic_cast<Flux::RangeExpression *>(&node)) {
    return eval(*range);
  }

  if (auto vector = dynamic_cast<Flux::VectorExpression *>(&node)) {
    return eval(*vector);
  }

  if (auto cell_selector =
          dynamic_cast<Flux::CellSelectorExpression *>(&node)) {
    return eval(*cell_selector);
  }

  if (auto binary = dynamic_cast<Flux::BinaryExpression *>(&node)) {
    return eval(*binary);
  }

  if (auto gen_func = dynamic_cast<Flux::GenFuncCallExpression *>(&node)) {
    return eval(*gen_func);
  }

  throw std::runtime_error("Unknown expression type in eval!");
}

Value Runtime::eval(Flux::LiteralExpression &node) {
  if (!_config.dim.has_value()) {
    throw std::runtime_error(
        "Dimension type must be specified before using LiteralExpression!");
  }

  if (_config.dim == Flux::DimType::Two) {
    float value = static_cast<float>(node.value);
    return Value2D([value](float i, float j, float t) { return value; });
  }

  if (_config.dim == Flux::DimType::Three) {
    float value = static_cast<float>(node.value);
    return Value3D(
        [value](float i, float j, float k, float t) { return value; });
  }

  throw std::runtime_error("Unknown dimension type in LiteralExpression!");
};

Value Runtime::eval(Flux::GenVariableExpression &node) {
  if (!_config.dim.has_value()) {
    throw std::runtime_error(
        "Dimension type must be specified before using GenVariableExpression!");
  }

  if (_config.dim == Flux::DimType::Two) {
    switch (node.var) {
    case Flux::GenVar::I:
      return Value2D([](float i, float j, float t) { return i; });
    case Flux::GenVar::J:
      return Value2D([](float i, float j, float t) { return j; });
    case Flux::GenVar::T:
      return Value2D([](float i, float j, float t) { return t; });
    default:
      throw std::runtime_error("Unsupported GenVar for 2D simulation!");
    }
  }

  if (_config.dim == Flux::DimType::Three) {
    switch (node.var) {
    case Flux::GenVar::I:
      return Value3D([](float i, float j, float k, float t) { return i; });
    case Flux::GenVar::J:
      return Value3D([](float i, float j, float k, float t) { return j; });
    case Flux::GenVar::K:
      return Value3D([](float i, float j, float k, float t) { return k; });
    case Flux::GenVar::T:
      return Value3D([](float i, float j, float k, float t) { return t; });
    default:
      throw std::runtime_error("Unsupported GenVar for 3D simulation!");
    }
  }

  throw std::runtime_error("Unknown dimension type in GenVariableExpression!");
};

Value Runtime::eval(Flux::RangeExpression &node) {
  if (!_config.dim.has_value()) {
    throw std::runtime_error(
        "Dimension type must be specified before using RangeExpression!");
  }

  auto lower_bound = eval(*node.start);
  auto upper_bound = eval(*node.end);

  if (_config.dim == Flux::DimType::Two) {
    auto lambda = [lower_bound = std::move(lower_bound),
                   upper_bound = std::move(upper_bound),
                   var = node.var](float i, float j, float t) {
      float lb = extract_float_2d(lower_bound, i, j, t);
      float ub = extract_float_2d(upper_bound, i, j, t);
      float value;
      switch (var) {
      case Flux::GenVar::I:
        value = i;
        break;
      case Flux::GenVar::J:
        value = j;
        break;
      case Flux::GenVar::K:
        throw std::runtime_error(
            "GenVar K is not valid in 2D RangeExpression!");
        break;
      case Flux::GenVar::T:
        value = t;
        break;
      }

      return (value >= lb && value <= ub) ? 1.0f : 0.0f;
    };
    return Value2D(std::move(lambda));
  }

  if (_config.dim == Flux::DimType::Three) {
    auto lambda = [lower_bound = std::move(lower_bound),
                   upper_bound = std::move(upper_bound),
                   var = node.var](float i, float j, float k, float t) {
      float lb = extract_float_3d(lower_bound, i, j, k, t);
      float ub = extract_float_3d(upper_bound, i, j, k, t);
      float value;
      switch (var) {
      case Flux::GenVar::I:
        value = i;
        break;
      case Flux::GenVar::J:
        value = j;
        break;
      case Flux::GenVar::K:
        value = k;
        break;
      case Flux::GenVar::T:
        value = t;
      }

      return (value >= lb && value <= ub) ? 1.0f : 0.0f;
    };
    return Value3D(std::move(lambda));
  }

  throw std::runtime_error("Unknown dimension type in RangeExpression!");
};

std::shared_ptr<std::vector<Value>>
Runtime::eval(Flux::VectorExpression &node) {
  std::vector<Value> components;
  for (auto &comp_expr : node.components) {
    components.push_back(std::get<Value>(eval(*comp_expr)));
  };

  return std::make_unique<std::vector<Value>>(std::move(components));
};

Membership Runtime::eval(Flux::CellSelectorExpression &node) {
  auto condition = eval(*node.condition);

  if (!_config.dim.has_value()) {
    throw std::runtime_error("Dimension type must be specified before using "
                             "CellSelectorExpression!");
  }

  if (_config.dim == Flux::DimType::Two) {
    return Membership2D(
        [condition = std::move(condition)](int i, int j, float t) {
          float cond_value;
          if (std::holds_alternative<Value>(condition)) {
            const auto &val = std::get<Value>(condition);
            cond_value = std::get<Value2D>(val)(static_cast<float>(i),
                                                static_cast<float>(j), t);
          } else {
            throw std::runtime_error(
                "Invalid condition type in CellSelectorExpression for 2D!");
          }

          return cond_value != 0.0f;
        });
  }

  if (_config.dim == Flux::DimType::Three) {
    return Membership3D([condition = std::move(condition)](int i, int j, int k,
                                                           float t) {
      float cond_value;
      if (std::holds_alternative<Value>(condition)) {
        const auto &val = std::get<Value>(condition);
        cond_value =
            std::get<Value3D>(val)(static_cast<float>(i), static_cast<float>(j),
                                   static_cast<float>(k), t);
      } else {
        throw std::runtime_error(
            "Invalid condition type in CellSelectorExpression for 3D!");
      }

      return cond_value != 0.0f;
    });
  }

  throw std::runtime_error("Unknown dimension type in CellSelectorExpression!");
};

Value Runtime::eval(Flux::BinaryExpression &node) {
  auto left = eval(*node.left);
  auto right = eval(*node.right);

  if (!_config.dim.has_value()) {
    throw std::runtime_error(
        "Dimension type must be specified before using BinaryExpression!");
  }

  if (_config.dim == Flux::DimType::Two) {
    auto lambda = [left = std::move(left), right = std::move(right),
                   op = node.op](float i, float j, float t) {
      float left_value = extract_float_2d(left, i, j, t);
      float right_value = extract_float_2d(right, i, j, t);

      switch (op) {
      case Flux::BinaryOp::Add:
        return left_value + right_value;
      case Flux::BinaryOp::Subtract:
        return left_value - right_value;
      case Flux::BinaryOp::Multiply:
        return left_value * right_value;
      case Flux::BinaryOp::Divide:
        return left_value / right_value;
      case Flux::BinaryOp::Power:
        return std::pow(left_value, right_value);
      case Flux::BinaryOp::And:
        return (left_value != 0.0f && right_value != 0.0f) ? 1.0f : 0.0f;
      case Flux::BinaryOp::Or:
        return (left_value != 0.0f || right_value != 0.0f) ? 1.0f : 0.0f;
      case Flux::BinaryOp::Equal:
        return (left_value == right_value) ? 1.0f : 0.0f;
      case Flux::BinaryOp::NotEqual:
        return (left_value != right_value) ? 1.0f : 0.0f;
      case Flux::BinaryOp::Less:
        return (left_value < right_value) ? 1.0f : 0.0f;
      case Flux::BinaryOp::LessEqual:
        return (left_value <= right_value) ? 1.0f : 0.0f;
      case Flux::BinaryOp::Greater:
        return (left_value > right_value) ? 1.0f : 0.0f;
      case Flux::BinaryOp::GreaterEqual:
        return (left_value >= right_value) ? 1.0f : 0.0f;
      default:
        throw std::runtime_error(
            "Unsupported BinaryOp in 2D BinaryExpression!");
      }
    };
    return Value2D(std::move(lambda));
  }

  if (_config.dim == Flux::DimType::Three) {
    auto lambda = [left = std::move(left), right = std::move(right),
                   op = node.op](float i, float j, float k, float t) {
      float left_value = extract_float_3d(left, i, j, k, t);
      float right_value = extract_float_3d(right, i, j, k, t);

      switch (op) {
      case Flux::BinaryOp::Add:
        return left_value + right_value;
      case Flux::BinaryOp::Subtract:
        return left_value - right_value;
      case Flux::BinaryOp::Multiply:
        return left_value * right_value;
      case Flux::BinaryOp::Divide:
        return left_value / right_value;
      case Flux::BinaryOp::Power:
        return std::pow(left_value, right_value);
      case Flux::BinaryOp::And:
        return (left_value != 0.0f && right_value != 0.0f) ? 1.0f : 0.0f;
      case Flux::BinaryOp::Or:
        return (left_value != 0.0f || right_value != 0.0f) ? 1.0f : 0.0f;
      case Flux::BinaryOp::Equal:
        return (left_value == right_value) ? 1.0f : 0.0f;
      case Flux::BinaryOp::NotEqual:
        return (left_value != right_value) ? 1.0f : 0.0f;
      case Flux::BinaryOp::Less:
        return (left_value < right_value) ? 1.0f : 0.0f;
      case Flux::BinaryOp::LessEqual:
        return (left_value <= right_value) ? 1.0f : 0.0f;
      case Flux::BinaryOp::Greater:
        return (left_value > right_value) ? 1.0f : 0.0f;
      case Flux::BinaryOp::GreaterEqual:
        return (left_value >= right_value) ? 1.0f : 0.0f;
      default:
        throw std::runtime_error(
            "Unsupported BinaryOp in 3D BinaryExpression!");
      }
    };
    return Value3D(std::move(lambda));
  }

  throw std::runtime_error("Unknown dimension type in BinaryExpression!");
}

Value Runtime::eval(Flux::GenFuncCallExpression &node) {
  std::vector<UnifiedValue> evaled_args;
  for (auto &arg_expr : node.arguments) {
    evaled_args.push_back(eval(*arg_expr));
  }

  if (!_config.dim.has_value()) {
    throw std::runtime_error(
        "Dimension type must be specified before using GenFuncCallExpression!");
  }

  if (_config.dim == Flux::DimType::Two) {
    auto lambda = [arguments = std::move(evaled_args),
                   func = node.func](float i, float j, float t) {
      std::vector<float> arg_values;
      for (const auto &arg : arguments) {
        float val = extract_float_2d(arg, i, j, t);
        arg_values.push_back(val);
      }

      switch (func) {
      case Flux::GenFunc::Sin:
        if (arg_values.size() != 1) {
          throw std::runtime_error(
              "Sin function requires exactly one argument!");
        }
        return std::sin(arg_values[0]);
      case Flux::GenFunc::Cos:
        if (arg_values.size() != 1) {
          throw std::runtime_error(
              "Cos function requires exactly one argument!");
        }
        return std::cos(arg_values[0]);
      case Flux::GenFunc::Tan:
        if (arg_values.size() != 1) {
          throw std::runtime_error(
              "Tan function requires exactly one argument!");
        }
        return std::tan(arg_values[0]);
      case Flux::GenFunc::Abs:
        if (arg_values.size() != 1) {
          throw std::runtime_error(
              "Abs function requires exactly one argument!");
        }
        return std::abs(arg_values[0]);
      case Flux::GenFunc::Sqrt:
        if (arg_values.size() != 1) {
          throw std::runtime_error(
              "Sqrt function requires exactly one argument!");
        }
        return std::sqrt(arg_values[0]);
      case Flux::GenFunc::Log:
        if (arg_values.size() != 1) {
          throw std::runtime_error(
              "Log function requires exactly one argument!");
        }
        return std::log(arg_values[0]);
      case Flux::GenFunc::Exp:
        if (arg_values.size() != 1) {
          throw std::runtime_error(
              "Exp function requires exactly one argument!");
        }
        return std::exp(arg_values[0]);
      case Flux::GenFunc::Min:
        if (arg_values.size() < 1) {
          throw std::runtime_error(
              "Min function requires at least one argument!");
        }
        return *std::min_element(arg_values.begin(), arg_values.end());
      case Flux::GenFunc::Max:
        if (arg_values.size() < 1) {
          throw std::runtime_error(
              "Max function requires at least one argument!");
        }
        return *std::max_element(arg_values.begin(), arg_values.end());
      default:
        throw std::runtime_error(
            "Unsupported GenFunc in 2D GenFuncCallExpression!");
      }
    };
    return Value2D(std::move(lambda));
  }

  if (_config.dim == Flux::DimType::Three) {
    auto lambda = [arguments = std::move(evaled_args),
                   func = node.func](float i, float j, float k, float t) {
      std::vector<float> arg_values;
      for (const auto &arg : arguments) {
        float val = extract_float_3d(arg, i, j, k, t);
        arg_values.push_back(val);
      }

      switch (func) {
      case Flux::GenFunc::Sin:
        if (arg_values.size() != 1) {
          throw std::runtime_error(
              "Sin function requires exactly one argument!");
        }
        return std::sin(arg_values[0]);
      case Flux::GenFunc::Cos:
        if (arg_values.size() != 1) {
          throw std::runtime_error(
              "Cos function requires exactly one argument!");
        }
        return std::cos(arg_values[0]);
      case Flux::GenFunc::Tan:
        if (arg_values.size() != 1) {
          throw std::runtime_error(
              "Tan function requires exactly one argument!");
        }
        return std::tan(arg_values[0]);
      case Flux::GenFunc::Abs:
        if (arg_values.size() != 1) {
          throw std::runtime_error(
              "Abs function requires exactly one argument!");
        }
        return std::abs(arg_values[0]);
      case Flux::GenFunc::Sqrt:
        if (arg_values.size() != 1) {
          throw std::runtime_error(
              "Sqrt function requires exactly one argument!");
        }
        return std::sqrt(arg_values[0]);
      case Flux::GenFunc::Log:
        if (arg_values.size() != 1) {
          throw std::runtime_error(
              "Log function requires exactly one argument!");
        }
        return std::log(arg_values[0]);
      case Flux::GenFunc::Exp:
        if (arg_values.size() != 1) {
          throw std::runtime_error(
              "Exp function requires exactly one argument!");
        }
        return std::exp(arg_values[0]);
      case Flux::GenFunc::Min:
        if (arg_values.size() < 1) {
          throw std::runtime_error(
              "Min function requires at least one argument!");
        }
        return *std::min_element(arg_values.begin(), arg_values.end());
      case Flux::GenFunc::Max:
        if (arg_values.size() < 1) {
          throw std::runtime_error(
              "Max function requires at least one argument!");
        }
        return *std::max_element(arg_values.begin(), arg_values.end());
      default:
        throw std::runtime_error(
            "Unsupported GenFunc in 2D GenFuncCallExpression!");
      }
    };

    return Value3D(std::move(lambda));
  }

  throw std::runtime_error("Unknown dimension type in GenFuncCallExpression!");
};

void Runtime::eval(Flux::Statement &node) {
  if (auto dim_stmt = dynamic_cast<Flux::DimStatement *>(&node)) {
    eval(*dim_stmt);
    return;
  }

  if (auto grid_stmt = dynamic_cast<Flux::GridStatement *>(&node)) {
    eval(*grid_stmt);
    return;
  }

  if (auto window_stmt = dynamic_cast<Flux::WindowStatement *>(&node)) {
    eval(*window_stmt);
    return;
  }

  if (auto set_stmt = dynamic_cast<Flux::SetStatement *>(&node)) {
    eval(*set_stmt);
    return;
  }

  if (auto density_stmt = dynamic_cast<Flux::DensityStatement *>(&node)) {
    eval(*density_stmt);
    return;
  }

  if (auto force_stmt = dynamic_cast<Flux::ForceStatement *>(&node)) {
    eval(*force_stmt);
    return;
  }

  if (auto solid_stmt = dynamic_cast<Flux::SolidStatement *>(&node)) {
    eval(*solid_stmt);
    return;
  }

  if (auto flow_stmt = dynamic_cast<Flux::FlowStatement *>(&node)) {
    eval(*flow_stmt);
    return;
  }

  if (auto flow_ratio_stmt = dynamic_cast<Flux::FlowRatioStatement *>(&node)) {
    eval(*flow_ratio_stmt);
    return;
  }

  throw std::runtime_error("Unknown statement type in eval!");
}

void Runtime::eval(Flux::DimStatement &node) {
  _config.dim = node.dim;
  if (_config.dim != Flux::DimType::Two &&
      _config.dim != Flux::DimType::Three) {
    throw std::runtime_error("Unsupported dimension type!");
  }
};

void Runtime::eval(Flux::GridStatement &node) {
  if (!_config.dim.has_value()) {
    throw std::runtime_error(
        "Dimension type must be specified before GridStatement!");
  }

  if (_config.dim == Flux::DimType::Two) {
    if (node.sizes.size() != 2) {
      throw std::runtime_error(
          "GridStatement size mismatch for 2D simulation!");
    }

    int nx = static_cast<int>(
        dynamic_cast<Flux::LiteralExpression *>(node.sizes[0].get())->value);
    int ny = static_cast<int>(
        dynamic_cast<Flux::LiteralExpression *>(node.sizes[1].get())->value);

    _config.simulation = std::make_unique<Simulation2D>(
        nx, ny, static_cast<float>(node.dx), static_cast<float>(node.dx));
  }

  if (_config.dim == Flux::DimType::Three) {
    if (node.sizes.size() != 3) {
      throw std::runtime_error(
          "GridStatement size mismatch for 3D simulation!");
    }

    int nx = static_cast<int>(
        dynamic_cast<Flux::LiteralExpression *>(node.sizes[0].get())->value);
    int ny = static_cast<int>(
        dynamic_cast<Flux::LiteralExpression *>(node.sizes[1].get())->value);
    int nz = static_cast<int>(
        dynamic_cast<Flux::LiteralExpression *>(node.sizes[2].get())->value);

    _config.simulation = std::make_unique<Simulation3D>(
        nx, ny, nz, static_cast<float>(node.dx), static_cast<float>(node.dx),
        static_cast<float>(node.dx));
  }
};

void Runtime::eval(Flux::WindowStatement &node) {
  if (!_config.dim.has_value()) {
    throw std::runtime_error(
        "Dimension type must be specified before WindowStatement!");
  }

  if (_config.dim == Flux::DimType::Two) {
    if (node.sizes.size() != 2) {
      throw std::runtime_error(
          "WindowStatement size mismatch for 2D simulation!");
    }

    float width =
        std::get<Value2D>(std::get<Value>(eval(*node.sizes[0])))(0, 0, 0);
    float height =
        std::get<Value2D>(std::get<Value>(eval(*node.sizes[1])))(0, 0, 0);

    auto &sim = std::get<std::unique_ptr<Simulation<2>>>(_config.simulation);

    _config.renderer = std::make_unique<Renderer2D>(
        dynamic_cast<Simulation2D &>(*sim).get_mac(), width, height);
  }

  if (_config.dim == Flux::DimType::Three) {
    if (node.sizes.size() != 2) {
      throw std::runtime_error(
          "WindowStatement size mismatch for 3D simulation!");
    }

    float width =
        std::get<Value3D>(std::get<Value>(eval(*node.sizes[0])))(0, 0, 0, 0);
    float height =
        std::get<Value3D>(std::get<Value>(eval(*node.sizes[1])))(0, 0, 0, 0);

    auto &sim = std::get<std::unique_ptr<Simulation<3>>>(_config.simulation);

    _config.renderer = std::make_unique<Renderer3D>(
        dynamic_cast<Simulation3D &>(*sim).get_mac(), width, height);
  }
};

void Runtime::eval(Flux::SetStatement &node) {
  auto membership = eval(*node.selector);
  _config.set_memberships[node.identifier] = std::move(membership);
  _config.set_order.push_back(node.identifier);
};

void Runtime::eval(Flux::DensityStatement &node) {
  auto value = std::get<Value>(eval(*node.value));
  _config.density_values[node.identifier] = std::move(value);
};

void Runtime::eval(Flux::ForceStatement &node) {
  auto value = eval(*node.value);
  _config.force_values[node.identifier] = std::move(value);
};

void Runtime::eval(Flux::SolidStatement &node) {
  auto membership = std::get<Value>(eval(*node.value));
  _config.solid_values[node.identifier] = std::move(membership);
};

void Runtime::eval(Flux::FlowStatement &node) {
  auto membership = std::get<Value>(eval(*node.value));
  _config.flow_values[node.identifier] = std::move(membership);
};

void Runtime::eval(Flux::FlowRatioStatement &node) {
  auto membership = std::get<Value>(eval(*node.value));
  _config.flow_ratio_values[node.identifier] = std::move(membership);
}
