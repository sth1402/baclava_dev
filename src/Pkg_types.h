#pragma once

namespace Pkg {
enum class Param { MU, RATE, SCALE, MEDIAN, UNK };
enum class Dist { WEIBULL, GAMMA, BERNOULLI, UNK };
enum class Model { CONSTANT, SIGMOID, LINEAR, EXP, STEP, UNK };
enum class EP { CENSORED, PRECLINICAL, CLINICAL };
enum class Step { MH, Gibbs, MH_with_Z };
enum class Components : uint8_t {
  UNK  = uint8_t(1) << 0,
  H    = uint8_t(1) << 1,
  P    = uint8_t(1) << 2,
  BETA = uint8_t(1) << 3,
  PSI  = uint8_t(1) << 4,
  Z    = uint8_t(1) << 5,
  TAU  = uint8_t(1) << 6
};
}