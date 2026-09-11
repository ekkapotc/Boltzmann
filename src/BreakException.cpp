#include "../inc/BreakException.hpp"

using namespace boltzmann;

const char* BreakException::what() const throw()
{
    return "Break Exception happened";
}
