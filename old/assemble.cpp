#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <tuple>
#include <functional>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/strong_components.hpp>

#include "khmer/khmer.hh"
#include "khmer/kmer_hash.hh"

struct Vertex
{
    khmer::HashIntoType kmer_hash;
    /* std::size_t index; */
    std::size_t component;
};

/* using VertexProperty = boost::property<boost::vertex_name_t, std::string>; */

using DiGraph = boost::adjacency_list<boost::vecS, boost::vecS,
                                      boost::bidirectionalS, Vertex>;

using vertex_dsc = DiGraph::vertex_descriptor;
using vertex_it = DiGraph::vertex_iterator;

using edge_dsc = DiGraph::edge_descriptor;
using edge_it = DiGraph::edge_iterator;

using Edge = std::pair<vertex_dsc, vertex_dsc>;

using ComponentDiGraph = 
            boost::filtered_graph<DiGraph, 
                                  std::function<bool(edge_dsc)>,
                                  std::function<bool(vertex_dsc)>>;

using KmerMap = std::unordered_map<khmer::HashIntoType, vertex_dsc>;

DiGraph create_graph(const std::string &infilename, khmer::WordLength k)
{
    std::ifstream in(infilename);

    if (!in) {
        std::cerr << "Error opening " << infilename << "\n";
        exit(1);
    }

    KmerMap kmerMap; // kmerMap[kmer] -> vertex in DG
    using hash_t = khmer::HashIntoType;
    using hash_pair = std::pair<hash_t, hash_t>;
    boost::unordered_set<hash_pair> hash_edges;
    boost::unordered_set<hash_t> kmer_hashes;

    std::cerr << "reading from file...\n\n";

    std::string read;
    // iterate over reads from input
    while (in >> read) {

        // ignore short reads
        if (read.size() <= k)
            continue;

        /* std::cerr << read << std::endl; */
        std::string kmer(read.begin(), read.begin()+k);
        hash_t hash_a = khmer::_hash_forward(kmer.c_str(), k), hash_b;

        kmer_hashes.emplace(hash_a);

        // iterate over kmers in read
        for (auto it = read.begin(); it+k != read.end(); ++it) {


            kmer = std::string(it+1, it+k+1);
            hash_b = khmer::_hash_forward(kmer.c_str(), k);

            /* std::cerr << kmer_a << " " << kmer_b << std::endl; */
            kmer_hashes.emplace(hash_b);
            hash_edges.emplace(hash_a, hash_b);

            /* kmer_a = kmer_b; */
            hash_a = hash_b;
        }
    }

    std::cerr << "adding vertices...\n\n";

    DiGraph DG;

    using v_prop = DiGraph::vertex_property_type;

    /* int i = 0; */
    for (auto &&hash : kmer_hashes) {
        v_prop v;
        v.kmer_hash = hash;
        /* v.index = i++; */

        kmerMap[hash] = add_vertex(v, DG);
    }

    std::cerr << "adding edges...\n\n";

    for (auto &&edge : hash_edges) {
        auto &&u = edge.first;
        auto &&v = edge.second;

        add_edge(kmerMap[u], kmerMap[v], DG);
    }

    return DG;
}

/* using mp = std::map<vertex_dsc, DiGraph::vertices_size_type>; */
using mp = std::vector<DiGraph::vertices_size_type>;
using bmp = boost::associative_property_map<mp>;
using val = mp::value_type;

std::tuple<val, mp>
weakly_connected_components(DiGraph &DG)
{
    /* std::vector<std::size_t> component(num_vertices(DG)); */

    std::vector<Edge> edges;

    vertex_dsc u_dsc, v_dsc;
    edge_it e_it, e_end;

    std::cerr << "saving original edges...\n\n";
    // iterate over edges of DG, add them to vector
    for (auto e_dsc : boost::make_iterator_range(boost::edges(DG))) {
        u_dsc = boost::source(e_dsc, DG);
        v_dsc = boost::target(e_dsc, DG);
        edges.emplace_back(u_dsc, v_dsc);
    }

    /* // iterate over edges of DG, add them to vector */
    /* for (std::tie(e_it, e_end) = boost::edges(DG); */
    /*         e_it != e_end; ++e_it) { */

    /*     edge_dsc &&e_dsc = *e_it; */

    /*     u_dsc = boost::source(e_dsc, DG); */
    /*     v_dsc = boost::target(e_dsc, DG); */
    /*     edges.emplace_back(u_dsc, v_dsc); */
    /* } */

    std::cerr << "adding reverse edges...\n\n";

    // for every edge in DG,
    // add its reverse to DG
    for (auto &&edge : edges) {
        vertex_dsc &u = edge.first;
        vertex_dsc &v = edge.second;
        boost::add_edge(v, u, DG);
    }

    std::cerr << "finding strongly connected components...\n\n";

    /* mp component; */
    /* bmp prop_map(component); */

    /* val num = boost::strong_components(DG, prop_map, vertex_index_map(get(&Vertex::index, DG))); */

    mp component(num_vertices(DG));

    val num = boost::strong_components(DG, component.data());

    std::cerr << "removing reverse edges...\n\n";

    std::reverse(edges.begin(), edges.end());

    // remove reverse edges from DG in reverse order
    for (auto &&edge : edges) {
        vertex_dsc &u = edge.first;
        vertex_dsc &v = edge.second;
        boost::remove_edge(v, u, DG);
    }

    return std::make_tuple(num, component);
}

std::vector<ComponentDiGraph>
create_component_digraphs(DiGraph &DG, khmer::WordLength k)
{
    val num;
    mp component;

    std::tie(num, component) = weakly_connected_components(DG);

    std::cerr << "total number of components: " << num << "\n\n";

    std::cerr << "finding component sizes...\n\n";

    std::vector<std::size_t> component_size(num);

    for (auto v_dsc : boost::make_iterator_range(vertices(DG))) {
        auto &&comp_num = component[v_dsc];
        DG[v_dsc].component = comp_num;
        ++component_size[comp_num];
    }

    std::cerr << "creating component filtered graphs...\n\n";

    std::vector<ComponentDiGraph> component_digraphs;

    for (std::size_t i = 0; i < num; ++i) {
        if (component_size[i] > 500u - k)
            component_digraphs.emplace_back(
                    DG,
                    [&DG, i](edge_dsc e) {
                        auto &&v = DG[source(e, DG)];
                        return v.component == i;
                    },
                    [&DG, i](vertex_dsc v) {
                        return DG[v].component == i;
                    });
    }

    return component_digraphs;
}

int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <reads> <k>\n";
        return 1;
    }

    std::string infilename(argv[1]);

    khmer::WordLength k = std::stoi(argv[2]);

    DiGraph DG = create_graph(infilename, k);

    std::cerr << "vertices: " << num_vertices(DG) << "\n";
    std::cerr << "edges:    "  << num_edges(DG) << "\n\n";

    std::cerr << "finding components...\n\n";

    std::vector<ComponentDiGraph> 
        component_digraphs = create_component_digraphs(DG, k);

    std::vector<std::string> assemblies;

    std::cerr << "iterating over component filtered graphs...\n\n";

    /* std::size_t total = component_digraphs.size(); */
    std::size_t semi_euler_cnt = 0;
    std::size_t non_semi_euler_cnt = 0;
    for (const auto &subgraph : component_digraphs) {

        /* std::cerr << total - ++i << " "; */

        std::size_t odd_cnt = 0;

        ComponentDiGraph::vertex_descriptor start, current, next;
        ComponentDiGraph::degree_size_type in_deg, out_deg;

        /* ComponentDiGraph::vertex_iterator v_it, v_end; */
        /* for (std::tie(v_it, v_end) = vertices(subgraph); */
        /*         v_it != v_end; ++v_it) { */
        /*     auto &&v_dsc = *v_it; */

        bool has_start = false, has_end = false;
        for (auto v_dsc : boost::make_iterator_range(vertices(subgraph))) {
            in_deg = in_degree(v_dsc, subgraph);
            out_deg = out_degree(v_dsc, subgraph);
            if (in_deg != out_deg) {
                ++odd_cnt;

                if (odd_cnt == 3)
                    break;

                if (out_deg == in_deg + 1) {
                    start = v_dsc;
                    has_start = true;
                }
                else if (out_deg + 1 == in_deg) {
                    has_end = true;
                }
            }
        }

        if (odd_cnt == 2 && has_start && has_end) {
            current = start;
            ++semi_euler_cnt;
        }
        else if (odd_cnt == 0) {
            current = *vertices(subgraph).first;
            ++semi_euler_cnt;
        }
        else {
            ++non_semi_euler_cnt;
            continue;
        }

        std::vector<ComponentDiGraph::vertex_descriptor> vertex_stack, path;

        /* ComponentDiGraph::edge_iterator out_edge_it; */

        vertex_stack.push_back(current);
        auto out_edge_dsc = *out_edges(current, subgraph).first;
        next = target(out_edge_dsc, subgraph);

        remove_edge(current, next, DG);
        current = next;

        while (!vertex_stack.empty()) {
            out_deg = out_degree(current, subgraph);

            if (out_deg == 0) {
                path.push_back(current);
                current = vertex_stack.back();
                vertex_stack.pop_back();
            }
            else {
                vertex_stack.push_back(current);
                out_edge_dsc = *out_edges(current, subgraph).first;
                next = target(out_edge_dsc, subgraph);

                remove_edge(current, next, DG);
                current = next;
            }
        }

        path.push_back(current);

        std::reverse(path.begin(), path.end());

        auto v = subgraph[path[0]];
        std::string assembly = khmer::_revhash(v.kmer_hash, k);

        for (std::size_t i = 1; i < path.size(); ++i) {
            v = subgraph[path[i]];

            auto tail = v.kmer_hash;

            assembly.push_back(khmer::_revhash(tail, 1)[0]);

            /* auto tail = v.kmer_hash << (k-1); */
            /* assembly.append(khmer::_revhash(tail, 1)); */
        }

        assemblies.push_back(assembly);

    }

    std::cerr << "Semi-eulerian subgraphs: " << semi_euler_cnt << std::endl
              << "Non-semi-eulerian:       " << non_semi_euler_cnt << std::endl;

    std::sort(assemblies.begin(), assemblies.end(),
              [](const std::string &str_a, const std::string &str_b)
                { return str_a.size() > str_b.size(); });

    std::size_t i = 0;
    for (auto &&assembly : assemblies) {
        std::cout << ">" << ++i << "\n";

        for (std::size_t i = 0; i < assembly.size(); i += 60) {
            std::cout << assembly.substr(i, 60) << "\n";
        }
    }

    return 0;
}
