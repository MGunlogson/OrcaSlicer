#ifndef slic3r_Magma_MagmaCell_hpp_
#define slic3r_Magma_MagmaCell_hpp_

#include <cstdint>
#include <functional>

namespace Slic3r {
namespace magma {

// Shape-agnostic cell id shared by all Magma lattices: up to three integer coordinates plus a
// `kind` tag for lattices with several cell types at one coordinate (tri-hex). What the fields
// mean is defined by each MagmaLattice implementation; for the triangle lattice (a, b, c) are
// the three triangle axes and `kind` is 0.
struct CellId {
    int     a = 0, b = 0, c = 0;
    uint8_t kind = 0;

    CellId() = default;
    CellId(int a_, int b_, int c_, uint8_t kind_ = 0)
        : a(a_), b(b_), c(c_), kind(kind_) {}

    bool operator==(const CellId &o) const {
        return a == o.a && b == o.b && c == o.c && kind == o.kind;
    }

    bool operator!=(const CellId &o) const {
        return !(*this == o);
    }

    // Lexicographic over (a, b, c, kind).
    bool operator<(const CellId &o) const {
        if (a != o.a) return a < o.a;
        if (b != o.b) return b < o.b;
        if (c != o.c) return c < o.c;
        return kind < o.kind;
    }
};

// Folds `kind` into the hash so cells that differ only by kind don't collide.
struct CellIdHash {
    size_t operator()(const CellId &c) const {
        size_t h = std::hash<int>()(c.a);
        h ^= std::hash<int>()(c.b) + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<int>()(c.c) + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<int>()(static_cast<int>(c.kind)) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

} // namespace magma
} // namespace Slic3r

#endif // slic3r_Magma_MagmaCell_hpp_
