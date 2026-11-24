#include "gradient.h"
#include "external/toml.hpp"

Gradient::Gradient(std::vector<double> c, std::vector<double> r,
                   std::vector<double> g, std::vector<double> b)
    : r(c, r), g(c, g), b(c, b) {}

png::rgb_pixel Gradient::operator()(double t) {
    return png::rgb_pixel(r(t), g(t), b(t));
}

std::vector<double> get_vector(const toml::table &tbl, std::string_view key) {
    const toml::array *arr = tbl[key].as_array();
    if (!arr) {
        throw std::runtime_error("Missing or invalid array for key: " +
                                 std::string(key));
    }
    std::vector<double> result;
    result.reserve(arr->size());
    for (size_t i = 0; i < arr->size(); ++i) {
        std::optional<double> val = (*arr)[i].value<double>();
        if (!val) {
            throw std::runtime_error("Non-numeric data found in '" +
                                     std::string(key) + "' at index " +
                                     std::to_string(i));
        }
        result.push_back(*val);
    }
    return result;
}

Gradient::Gradient(std::string filename) {
    toml::table tbl;
    try {
        tbl = toml::parse_file(filename);
    } catch (const toml::parse_error &err) {
        throw std::runtime_error(
            std::string("Bad TOML file passed to Gradient::Gradient.") +
            err.what());
    }
    const toml::table *gradient = tbl["gradient"].as_table();
    if (!gradient) {
        throw std::runtime_error("Missing or invalid [gradient] table.");
    }
    std::vector<double> c, r, g, b;
    try {
        c = get_vector(*gradient, "c");
        r = get_vector(*gradient, "r");
        g = get_vector(*gradient, "g");
        b = get_vector(*gradient, "b");
    } catch (const std::runtime_error &err) {
        throw std::runtime_error(
            std::string("Gradient must contain c, r, g, and b arrays:") +
            err.what());
    }
    if (c.empty() || r.size() != c.size() || g.size() != c.size() ||
        b.size() != c.size()) {
        throw std::runtime_error(
            "Gradient arrays must be non-empty and of equal length.");
    }
    this->r = Pchip(c, r);
    this->g = Pchip(c, g);
    this->b = Pchip(c, b);
}
