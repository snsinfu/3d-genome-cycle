#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <vector>

#include "io.hpp"


static std::unique_ptr<std::FILE, int(*)(std::FILE*)>
fopen_raii(char const* filename, char const* mode)
{
    return {std::fopen(filename, mode), std::fclose};
}


std::string
load_text(std::string const& filename)
{
    constexpr std::size_t buffer_size = 64 * 1024;

    auto file = fopen_raii(filename.c_str(), "rb");
    if (!file) {
        throw std::runtime_error("failed to open file " + filename);
    }

    std::string content;
    std::vector<char> buffer(buffer_size);

    for (;;) {
        auto const n_read = std::fread(buffer.data(), 1, buffer.size(), file.get());

        if (n_read > 0) {
            content.append(buffer.data(), std::size_t(n_read));
        }

        if (n_read < buffer.size()) {
            if (std::ferror(file.get())) {
                throw std::runtime_error("failed to read from file " + filename);
            }
            break;
        }
    }

    return content;
}
