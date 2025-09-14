// EdgesCGALDelaunay3D.cpp (clean warnings)
// Delaunay triangulation 3D edges extractor using CGAL, fast I/O and optional TBB parallelism.
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_data_structure_3.h>

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <limits>
#include <utility>
#include <algorithm>
#include <thread>

#ifdef CGAL_LINKED_WITH_TBB
  #include <tbb/global_control.h>
#endif

struct BBox {
    double xmin =  std::numeric_limits<double>::infinity();
    double ymin =  std::numeric_limits<double>::infinity();
    double zmin =  std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();
    double ymax = -std::numeric_limits<double>::infinity();
    double zmax = -std::numeric_limits<double>::infinity();
    void update(double x, double y, double z) {
        xmin = std::min(xmin, x); ymin = std::min(ymin, y); zmin = std::min(zmin, z);
        xmax = std::max(xmax, x); ymax = std::max(ymax, y); zmax = std::max(zmax, z);
    }
    bool valid() const { return xmin <= xmax && ymin <= ymax && zmin <= zmax; }
};

static bool load_npy_points(const char* path,
    std::vector<std::pair<CGAL::Exact_predicates_inexact_constructions_kernel::Point_3, uint32_t>>& pts,
    BBox& bb)
{
    using Point = CGAL::Exact_predicates_inexact_constructions_kernel::Point_3;
    std::FILE* f = std::fopen(path, "rb");
    if (!f) { std::perror("fopen(input)"); return false; }

    char magic[6];
    if (std::fread(magic,1,6,f) != 6 || std::memcmp(magic, "\x93NUMPY",6) != 0) {
        std::fprintf(stderr, "Invalid NPY magic string.\n");
        std::fclose(f); return false;
    }
    unsigned char ver[2];
    if (std::fread(ver,1,2,f) != 2) { std::fclose(f); return false; }
    uint32_t header_len = 0;
    if (ver[0] == 1) {
        uint16_t hl; if (std::fread(&hl,2,1,f) != 1) { std::fclose(f); return false; }
        header_len = hl;
    } else if (ver[0] == 2) {
        uint32_t hl; if (std::fread(&hl,4,1,f) != 1) { std::fclose(f); return false; }
        header_len = hl;
    } else {
        std::fprintf(stderr, "Unsupported NPY version %u.%u\n", ver[0], ver[1]);
        std::fclose(f); return false;
    }
    std::string header(header_len, '\0');
    if (std::fread(header.data(),1,header_len,f) != header_len) { std::fclose(f); return false; }
    if (header.find("True") != std::string::npos && header.find("fortran_order") != std::string::npos) {
        std::fprintf(stderr, "Fortran-ordered arrays not supported.\n");
        std::fclose(f); return false;
    }
    // descr
    size_t pos = header.find("'descr':");
    if (pos == std::string::npos) { std::fclose(f); return false; }
    size_t start = header.find("'", pos+8);
    size_t end = header.find("'", start+1);
    std::string descr = header.substr(start+1, end-start-1);
    bool is_f8 = (descr == "<f8");
    bool is_f4 = (descr == "<f4");
    if (!is_f8 && !is_f4) {
        std::fprintf(stderr, "Only little-endian float32/float64 supported.\n");
        std::fclose(f); return false;
    }
    // shape
    pos = header.find("'shape':");
    if (pos == std::string::npos) { std::fclose(f); return false; }
    start = header.find('(', pos);
    end = header.find(')', start);
    std::string shape = header.substr(start+1, end-start-1);
    const char* p = shape.c_str();
    char* q;
    size_t npts = std::strtoull(p, &q, 10);
    p = q;
    while (*p == ' ' || *p == ',') ++p;
    size_t dim = std::strtoull(p, &q, 10);
    if (dim != 3) {
        std::fprintf(stderr, "Input array must be of shape (N,3).\n");
        std::fclose(f); return false;
    }

    size_t total = npts * 3;
    pts.reserve(npts);
    if (is_f8) {
        std::vector<double> buf(total);
        if (std::fread(buf.data(), sizeof(double), total, f) != total) { std::fclose(f); return false; }
        for (size_t i = 0; i < npts; ++i) {
            double x = buf[3*i], y = buf[3*i+1], z = buf[3*i+2];
            if (std::isfinite(x) && std::isfinite(y) && std::isfinite(z)) {
                bb.update(x,y,z);
                pts.emplace_back(Point(x,y,z), (uint32_t)pts.size());
            }
        }
    } else {
        std::vector<float> buf(total);
        if (std::fread(buf.data(), sizeof(float), total, f) != total) { std::fclose(f); return false; }
        for (size_t i = 0; i < npts; ++i) {
            double x = buf[3*i], y = buf[3*i+1], z = buf[3*i+2];
            if (std::isfinite(x) && std::isfinite(y) && std::isfinite(z)) {
                bb.update(x,y,z);
                pts.emplace_back(Point(x,y,z), (uint32_t)pts.size());
            }
        }
    }
    std::fclose(f);
    return true;
}

static void write_npy_header(std::FILE* f, size_t rows) {
    std::string header = "{'descr': '<u4', 'fortran_order': False, 'shape': (" +
        std::to_string(rows) + ", 2), }";
    size_t header_len = header.size() + 1; // for trailing newline
    size_t pad = 16 - ((10 + header_len) % 16);
    header.append(pad, ' ');
    header.push_back('\n');
    uint16_t hlen = (uint16_t)header.size();
    std::fwrite("\x93NUMPY",1,6,f);
    unsigned char ver[2] = {1,0};
    std::fwrite(ver,1,2,f);
    std::fwrite(&hlen,2,1,f);
    std::fwrite(header.data(),1,header.size(),f);
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::fprintf(stderr,
            "Usage: %s input.npy output.npy\n"
            "input: (N,3) float32/float64 array\n", argv[0]);
        return 1;
    }
    const char* in_path = argv[1];
    const char* out_path = argv[2];

    using K   = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Vb0 = CGAL::Triangulation_vertex_base_3<K>;
    using Vb  = CGAL::Triangulation_vertex_base_with_info_3<uint32_t, K, Vb0>;
    using Cb  = CGAL::Delaunay_triangulation_cell_base_3<K>;
    #ifdef CGAL_LINKED_WITH_TBB
      using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb, CGAL::Parallel_tag>;
    #else
      using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb>;
    #endif
    using Dt  = CGAL::Delaunay_triangulation_3<K, Tds>;
    using Point = K::Point_3;

    #ifdef CGAL_LINKED_WITH_TBB
    {
        int threads = (int)std::thread::hardware_concurrency();
        if (const char* env = std::getenv("CGAL_NTHREADS")) {
            int req = std::atoi(env);
            if (req > 0) threads = req;
        }
        static tbb::global_control gc(tbb::global_control::max_allowed_parallelism, threads);
        std::fprintf(stderr, "[info] TBB enabled, using up to %d threads\n", threads);
    }
    #endif

    std::vector<std::pair<Point, uint32_t>> pts;
    BBox bb;
    if (!load_npy_points(in_path, pts, bb)) return 2;

    if (pts.size() < 2) {
        std::fprintf(stderr, "Input has < 2 valid points. Nothing to do.\n");
        std::FILE* fout = std::fopen(out_path, "wb");
        if (fout) std::fclose(fout);
        return 0;
    }
    std::fprintf(stderr, "[info] Loaded %zu points\n", pts.size());

    // Triangulation
    #ifdef CGAL_LINKED_WITH_TBB
    std::unique_ptr<Dt::Lock_data_structure> lock_holder;
    Dt::Lock_data_structure* lock_ptr = nullptr;
    if (bb.valid()) {
        CGAL::Bbox_3 box(bb.xmin, bb.ymin, bb.zmin, bb.xmax, bb.ymax, bb.zmax);
        lock_holder.reset(new Dt::Lock_data_structure(box, 64));
        lock_ptr = lock_holder.get();
    }
    Dt dt(Dt::Geom_traits(), lock_ptr);
    #else
    Dt dt;
    #endif

    dt.insert(pts.begin(), pts.end());
    if (!dt.is_valid(false))
        std::fprintf(stderr, "[warning] CGAL triangulation reports invalid structure.\n");

    // Count edges
    size_t edge_count = 0;
    for (auto eit = dt.finite_edges_begin(); eit != dt.finite_edges_end(); ++eit) {
        auto c = eit->first;
        int i = eit->second, j = eit->third;
        auto vi = c->vertex(i);
        auto vj = c->vertex(j);
        if (dt.is_infinite(vi) || dt.is_infinite(vj)) continue;
        uint32_t a = vi->info(), b = vj->info();
        if (a == b) continue;
        if (a > b) std::swap(a,b);
        ++edge_count;
    }

    std::FILE* fout = std::fopen(out_path, "wb");
    if (!fout) { std::perror("fopen(output)"); return 3; }
    write_npy_header(fout, edge_count);
    for (auto eit = dt.finite_edges_begin(); eit != dt.finite_edges_end(); ++eit) {
        auto c = eit->first;
        int i = eit->second, j = eit->third;
        auto vi = c->vertex(i);
        auto vj = c->vertex(j);
        if (dt.is_infinite(vi) || dt.is_infinite(vj)) continue;
        uint32_t a = vi->info(), b = vj->info();
        if (a == b) continue;
        if (a > b) std::swap(a,b);
        uint32_t ab[2] = {a,b};
        std::fwrite(ab, sizeof(uint32_t), 2, fout);
    }
    std::fclose(fout);
    std::fprintf(stderr, "[info] Wrote %zu edges to %s\n", edge_count, out_path);
    return 0;
}
