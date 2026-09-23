#ifndef NEGATIVERATEERROR_H
#define NEGATIVERATEERROR_H

#include <exception>
#include <string>

class NegativeRateError : public std::exception {
public:

    NegativeRateError(const std::string& message) : message_(message) {}
  
    const char* what() const noexcept override {
        return message_.c_str();
    }

private:
    std::string message_;
};

#endif
