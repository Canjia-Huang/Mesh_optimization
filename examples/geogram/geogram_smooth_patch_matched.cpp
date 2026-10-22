
#include <Mesh_smoothing_3/Mesh_smoothing_3.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/mesh/mesh_geometry.h>

#include <memory>
#include <vector>

class Mesh_wrapper : public Mesh_smoothing_3::helper_structures::Mixed_element_mesh<GEO::index_t, GEO::index_t, GEO::vec3, GEO::index_range> {
public:
    std::size_t nb_vertices() const override { return mesh.vertices.nb(); }

    Point_3 vertex_coordinates(Vertex_descriptor vertex) const override {
        return mesh.vertices.point(vertex);
    }
    void set_new_vertex_coordinates(GEO::index_t vertex, Point_3 coord) override { mesh.vertices.point(vertex) = coord; }

    GEO::index_range input_cell_range() const override { return GEO::index_range{mesh.cells.begin(), mesh.cells.end()}; }

    Vertex_descriptor get_cell_vertex(GEO::index_t cell, unsigned local_Vertex_descriptor) const override {
        return mesh.cells.vertex(cell, local_Vertex_descriptor);
    };

    GEO::vec3 get_ref_vertex_coordinates(GEO::index_t vertex) const override {
        assert(ref_mesh != nullptr);
        return ref_mesh->vertices.point(vertex);
    }

    Shape const * get_element_shape(GEO::index_t cell) const override {
        switch (mesh.cells.type(cell)) {
            case GEO::MESH_TET:
                return &tet_ref;
                break;
            case GEO::MESH_HEX:
                return &hex_ref;
                break;
            case GEO::MESH_PYRAMID:
                return &py_ref;
                break;
            case GEO::MESH_PRISM:
                return &we_ref;
                break;
            default:
                return nullptr;
                break;
        }
    }

public:
    Mesh_wrapper(GEO::Mesh &mesh_, GEO::Mesh const * ref_mesh_)
    : mesh(mesh_)
    , ref_mesh(ref_mesh_)
    {
        this->has_reference_mesh = (ref_mesh != nullptr);
        we_ref.inverse = true;
        this->assemble(); // CRITICAL
    }

    void set_orientation(bool inv_tet, bool inv_hex, bool inv_pyr, bool inv_wed) {
        tet_ref.inverse = inv_tet;
        hex_ref.inverse = inv_hex;
        py_ref.inverse  = inv_pyr;
        we_ref.inverse  = !inv_wed;
    }

    GEO::Mesh &mesh;
    GEO::Mesh const * ref_mesh;

    Mesh_smoothing_3::Shapes::VTK_TETRAHEDRON<GEO::vec3> tet_ref;
    Mesh_smoothing_3::Shapes::GEOGRAM_HEXAHEDRON<GEO::vec3> hex_ref;
    Mesh_smoothing_3::Shapes::VTK_PYRAMID<GEO::vec3> py_ref;
    Mesh_smoothing_3::Shapes::VTK_WEDGE<GEO::vec3> we_ref;
};

class Boundary_wrapper {
public:
    using Face_descriptor = GEO::index_t;
    using Normal_3 = GEO::vec3;
    using Surface_patch_index = unsigned;
    unsigned nb_faces() const { return mesh.facets.nb(); }

    GEO::index_range face_range() const { return GEO::index_range{mesh.facets.begin(), mesh.facets.end()}; }

    unsigned patch_id(Face_descriptor f) const {
        return facet_region.is_bound() ? static_cast<unsigned>(facet_region[f]) : 0;
    }
    unsigned nb_face_vertices(Face_descriptor face) const { return mesh.facets.nb_vertices(face); }
    auto face_vertices(Face_descriptor face) const { return mesh.facets.vertices(face); }
public:
    GEO::Mesh const &mesh;
    GEO::Attribute<GEO::index_t> facet_region;

    Boundary_wrapper(GEO::Mesh const &m) : mesh(m) {
        if (const_cast<GEO::Mesh&>(mesh).facets.attributes().is_defined("region")) {
            facet_region.bind(const_cast<GEO::Mesh&>(mesh).facets.attributes(), "region");
        } else {
            facet_region.bind_if_is_defined(const_cast<GEO::Mesh&>(mesh).facets.attributes(), "facet_region");
        }
    }
};

class PolyLine_wrapper {
public:
    using Edge_descriptor = int;
    using Curve_index = unsigned;
    unsigned nb_edges() const { return mesh.edges.nb(); }
    GEO::index_range edge_range() const {
        return GEO::index_range{mesh.edges.begin(), mesh.edges.end()};
    }
    unsigned curve_id(Edge_descriptor e) const {
        return edge_region.is_bound() ? static_cast<unsigned>(edge_region[e]) : 0;
    }
    int edge_vertex(Edge_descriptor edge, unsigned i) const { return mesh.edges.vertex(edge, i); }
public:
    GEO::Mesh const &mesh;
    GEO::Attribute<GEO::index_t> edge_region;

    PolyLine_wrapper(GEO::Mesh const &m) : mesh(m) {
        if (const_cast<GEO::Mesh&>(mesh).edges.attributes().is_defined("region")) {
            edge_region.bind(const_cast<GEO::Mesh&>(mesh).edges.attributes(), "region");
        } else {
            edge_region.bind_if_is_defined(const_cast<GEO::Mesh&>(mesh).edges.attributes(), "edge_region");
        }
    }
};

struct PatchData {
    GEO::Mesh mesh;
    std::unique_ptr<GEO::MeshFacetsAABB> aabb;
    std::vector<GEO::vec3> normal;
};

struct CurveData {
    GEO::Mesh mesh;
    std::unique_ptr<GEO::MeshFacetsAABB> aabb;
    std::vector<GEO::vec3> direction;
};

int main(int argc, char** argv) {
    GEO::initialize();
    const std::string filename = (argc > 1) ? argv[1] : "../data/fandisk_kenshi_hexmesh.mesh";
    const std::string surface_filename = (argc > 2) ? argv[2] : "";
    const std::string curves_filename = (argc > 3) ? argv[3] : "";
    const std::string point_targets_filename = (argc > 4) ? argv[4] : "";
    const std::string output_filename = (argc > 5) ? argv[5] : "output.mesh";

    GEO::MeshIOFlags ioflags;
    ioflags.set_attribute(GEO::MESH_FACET_REGION);
    ioflags.set_attribute(GEO::MESH_EDGE_REGION);

    GEO::Mesh mesh;
    if(!GEO::mesh_load(filename, mesh, ioflags)) {
        std::cerr << "Error loading mesh: " << filename << std::endl;
        return EXIT_FAILURE;
    }

    GEO::Mesh reference;
    reference.copy(mesh);

    GEO::Mesh surface;
    if (!surface_filename.empty()) {
        if(!GEO::mesh_load(surface_filename, surface, ioflags)) {
            std::cerr << "Error loading surface mesh: " << surface_filename << std::endl;
            return EXIT_FAILURE;
        }
    }
    else {
        surface.copy(mesh);
    }
    surface.edges.clear();
    surface.cells.clear();
    surface.vertices.remove_isolated();
    surface.facets.triangulate();

    GEO::Attribute<GEO::index_t> surf_facet_region;
    if (surface.facets.attributes().is_defined("region")) {
        surf_facet_region.bind(surface.facets.attributes(), "region");
    } else {
        surf_facet_region.bind_if_is_defined(surface.facets.attributes(), "facet_region");
    }

    // Partition surface facets into patch_trees by patch_id
    unsigned max_patch_id = 0;
    if (surf_facet_region.is_bound()) {
        for (auto f : surface.facets) {
            unsigned pid = static_cast<unsigned>(surf_facet_region[f]);
            if (pid > max_patch_id) max_patch_id = pid;
        }
    }

    std::vector<std::vector<GEO::index_t>> patch_facets(max_patch_id + 1);
    for (auto f : surface.facets) {
        unsigned pid = surf_facet_region.is_bound() ? static_cast<unsigned>(surf_facet_region[f]) : 0;
        patch_facets[pid].push_back(f);
    }

    std::vector<PatchData> patch_trees(max_patch_id + 1);
    std::vector<GEO::index_t> g2l(surface.vertices.nb(), GEO::NO_INDEX);

    for (unsigned pid = 0; pid <= max_patch_id; ++pid) {
        if (patch_facets[pid].empty()) continue;
        PatchData& pdata = patch_trees[pid];
        for (auto f : patch_facets[pid]) {
            GEO::index_t lvs[3];
            for (GEO::index_t i = 0; i < 3; ++i) {
                GEO::index_t gv = surface.facets.vertex(f, i);
                if (g2l[gv] == GEO::NO_INDEX) {
                    g2l[gv] = pdata.mesh.vertices.create_vertex();
                    pdata.mesh.vertices.point(g2l[gv]) = surface.vertices.point(gv);
                }
                lvs[i] = g2l[gv];
            }
            pdata.mesh.facets.create_triangle(lvs[0], lvs[1], lvs[2]);
        }
        for (auto f : patch_facets[pid]) {
            for (GEO::index_t i = 0; i < 3; ++i) {
                g2l[surface.facets.vertex(f, i)] = GEO::NO_INDEX;
            }
        }

        pdata.aabb = std::make_unique<GEO::MeshFacetsAABB>();
        pdata.aabb->initialize(pdata.mesh);
        pdata.normal.resize(pdata.mesh.facets.nb());
        for (auto f : pdata.mesh.facets) {
            pdata.normal[f] = GEO::Geom::mesh_facet_normal(pdata.mesh, f);
            if (pdata.normal[f].length() > 1e-10) pdata.normal[f] /= pdata.normal[f].length();
            else pdata.normal[f] = GEO::vec3(0., 0., 1.);
        }
    }

    const PatchData* default_patch = nullptr;
    for (auto& pdata : patch_trees) {
        if (pdata.aabb) {
            if (!default_patch || &pdata == &patch_trees[0]) default_patch = &pdata;
        }
    }

    auto query = [&](GEO::vec3 const &coord, unsigned surface_id, double radius) -> std::tuple<GEO::vec3, GEO::vec3, double> {
        const PatchData& p = (surface_id < patch_trees.size() && patch_trees[surface_id].aabb)
                           ? patch_trees[surface_id]
                           : *default_patch;
        GEO::vec3 res;
        double sqr_dist;
        GEO::index_t facet = p.aabb->nearest_facet(coord, res, sqr_dist);
        return {res, p.normal[facet], 1.};
    };


    GEO::Mesh curves;
    if (!curves_filename.empty()) {
        if(!GEO::mesh_load(curves_filename, curves, ioflags)) {
            std::cerr << "Error loading curves mesh: " << curves_filename << std::endl;
            return EXIT_FAILURE;
        }
    }
    else {
        curves.copy(mesh);
    }
    curves.facets.clear();
    curves.cells.clear();
    curves.vertices.remove_isolated();

    GEO::Attribute<GEO::index_t> curves_edge_region;
    if (curves.edges.attributes().is_defined("region")) {
        curves_edge_region.bind(curves.edges.attributes(), "region");
    } else {
        curves_edge_region.bind_if_is_defined(curves.edges.attributes(), "edge_region");
    }

    // Partition curve edges into curve_trees by curve_id
    unsigned max_curve_id = 0;
    if (curves_edge_region.is_bound()) {
        for (auto e : curves.edges) {
            unsigned cid = static_cast<unsigned>(curves_edge_region[e]);
            if (cid > max_curve_id) max_curve_id = cid;
        }
    }

    std::vector<std::vector<GEO::index_t>> curve_edges(max_curve_id + 1);
    for (auto e : curves.edges) {
        unsigned cid = curves_edge_region.is_bound() ? static_cast<unsigned>(curves_edge_region[e]) : 0;
        curve_edges[cid].push_back(e);
    }

    std::vector<CurveData> curve_trees(max_curve_id + 1);
    std::vector<GEO::index_t> c_g2l(curves.vertices.nb(), GEO::NO_INDEX);

    for (unsigned cid = 0; cid <= max_curve_id; ++cid) {
        if (curve_edges[cid].empty()) continue;
        CurveData& cdata = curve_trees[cid];
        cdata.direction.reserve(curve_edges[cid].size());
        for (auto e : curve_edges[cid]) {
            GEO::index_t gv0 = curves.edges.vertex(e, 0);
            GEO::index_t gv1 = curves.edges.vertex(e, 1);
            if (c_g2l[gv0] == GEO::NO_INDEX) {
                c_g2l[gv0] = cdata.mesh.vertices.create_vertex();
                cdata.mesh.vertices.point(c_g2l[gv0]) = curves.vertices.point(gv0);
            }
            if (c_g2l[gv1] == GEO::NO_INDEX) {
                c_g2l[gv1] = cdata.mesh.vertices.create_vertex();
                cdata.mesh.vertices.point(c_g2l[gv1]) = curves.vertices.point(gv1);
            }
            GEO::index_t lv0 = c_g2l[gv0];
            GEO::index_t lv1 = c_g2l[gv1];
            cdata.mesh.facets.create_triangle(lv0, lv1, lv1);
            GEO::vec3 dir = cdata.mesh.vertices.point(lv1) - cdata.mesh.vertices.point(lv0);
            if (dir.length() > 1e-10) dir /= dir.length();
            else dir = GEO::vec3(0., 0., 1.);
            cdata.direction.push_back(dir);
        }
        for (auto e : curve_edges[cid]) {
            c_g2l[curves.edges.vertex(e, 0)] = GEO::NO_INDEX;
            c_g2l[curves.edges.vertex(e, 1)] = GEO::NO_INDEX;
        }

        cdata.aabb = std::make_unique<GEO::MeshFacetsAABB>();
        cdata.aabb->initialize(cdata.mesh);
    }

    const CurveData* default_curve = nullptr;
    for (auto& cdata : curve_trees) {
        if (cdata.aabb) {
            if (!default_curve || &cdata == &curve_trees[0]) default_curve = &cdata;
        }
    }

    auto curve_query = [&](GEO::vec3 const &coord, unsigned edge_id, double radius) -> std::tuple<GEO::vec3, GEO::vec3, double> {
        if (default_curve) {
            const CurveData& c = (edge_id < curve_trees.size() && curve_trees[edge_id].aabb)
                               ? curve_trees[edge_id]
                               : *default_curve;
            GEO::vec3 res;
            double sqr_dist;
            GEO::index_t edge = c.aabb->nearest_facet(coord, res, sqr_dist);
            return {res, c.direction[edge], 1.};
        }
        return {coord, GEO::vec3(0., 0., 1.), 0.};
    };

    std::vector<std::tuple<GEO::index_t, GEO::vec3, double>> point_targets;

    if (!point_targets_filename.empty()) {
        std::ifstream infile(point_targets_filename);
        if (!infile) {
            std::cerr << "Error opening point targets file: " << point_targets_filename << std::endl;
            return EXIT_FAILURE;
        }
        GEO::index_t vertex_id;
        double x, y, z;
        while (infile >> vertex_id >> x >> y >> z) {
            point_targets.emplace_back(vertex_id, GEO::vec3(x, y, z), 1.);
        }
    }

    Mesh_wrapper mesh_wrapper(mesh, nullptr); // replacing nullptr by &reference will use it as a reference;
    mesh_wrapper.set_orientation(false, true, false, true); // beware of the orientation of your input elements!

    Boundary_wrapper surface_wrapper{mesh};
    PolyLine_wrapper curve_wrapper{mesh};

    Mesh_smoothing_3::Mesh_smoother smoother(mesh_wrapper, surface_wrapper, curve_wrapper);

    smoother.set_boundary_query(query);
    smoother.set_curves_query(curve_query);
    smoother.set_vertex_target_positions(point_targets);

    smoother.set_boundary_weight(Mesh_smoothing_3::Parameters::STRONG);
    smoother.set_verbose();
    smoother.set_max_number_of_iteration(100);
    smoother.run();

    GEO::mesh_save(mesh, output_filename);

    return EXIT_SUCCESS;
}

