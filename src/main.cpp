// Command-line front end for the Sturm root bracketing library.
// All console I/O lives here; the library itself is free of I/O.

#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "sturm/sturm.hpp"

namespace {

// Read a value of type T, re-prompting until the input is valid.
template <typename T>
T prompt(const std::string& message) {
    T value{};
    while (true) {
        std::cout << message;
        if (std::cin >> value) return value;
        if (std::cin.eof()) throw std::runtime_error("unexpected end of input");
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout << "  Invalid input, please try again.\n";
    }
}

void show_sequence(const std::vector<std::vector<double>>& seq) {
    std::cout << "\nSturm sequence:\n" << std::fixed;
    for (const auto& row : seq) {
        for (double c : row) std::cout << c << ' ';
        std::cout << '\n';
    }
}

void show_roots(const std::vector<double>& roots) {
    if (roots.empty()) {
        std::cout << "\nNo distinct real roots.\n";
    } else if (roots.size() == 1) {
        std::cout << "\n1 distinct real root at x = " << roots.front() << '\n';
    } else {
        std::cout << '\n' << roots.size() << " distinct real roots for x:";
        for (double r : roots) std::cout << ' ' << r;
        std::cout << '\n';
    }
}

}  // namespace

int main() {
    try {
        const int order = prompt<int>("Please enter order of polynomial: ");
        if (order < 0) {
            std::cerr << "Order must be non-negative.\n";
            return 1;
        }

        std::cout << '\n';
        std::vector<double> descending;  // coefficient of x^order first
        descending.reserve(static_cast<std::size_t>(order) + 1);
        for (int i = order; i >= 0; --i)
            descending.push_back(
                prompt<double>("Please enter coefficient number " + std::to_string(i) + ": "));

        const sturm::Sturm poly = sturm::Sturm::from_descending(descending);
        show_sequence(poly.sequence());
        show_roots(poly.real_roots());
        return 0;
    } catch (const sturm::invalid_polynomial& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}
