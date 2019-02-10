#define _CRT_SECURE_NO_WARNINGS 1

#include <fstream>
#include <utility>
#include <cstddef>
#include <stack>
#include <cctype>
#include <vector>
#include <optional>
#include <cassert>
#include <string_view>
#include <map>
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>

// vmf parsing

class vmf_group {
public:
	std::multimap<std::string, vmf_group> children;
	std::map<std::string, std::string> kvs;
};

void ignore_spaces(std::string_view data, std::size_t& offset) {
	while (offset < data.size() && std::isspace(data[offset])) {
		++offset;
	}
}

struct Token {
	enum class Kind {
		Identifier,
		String,
		LeftBrace,
		RightBrace,
	} kind;
	std::string value;
};

std::optional<Token> next_token(const std::string& data, std::size_t& offset) {
	ignore_spaces(data, offset);
	if (offset >= data.size()) { return std::nullopt; }
	std::size_t end = 0;
	std::string str;
	switch (data[offset]) {
	case '{':
		++offset;
		return Token{ Token::Kind::LeftBrace, "{" };
	case '}':
		++offset;
		return Token{ Token::Kind::RightBrace, "}" };
	case '"':
		end = data.find('"', offset + 1);
		assert(end != std::string::npos);
		str = data.substr(offset + 1, end - offset - 1);
		offset = end + 1;
		return Token{ Token::Kind::String, str };
	default:
		assert(std::isalnum(data[offset]) || data[offset] == '_');
		while (offset <= data.size() && (std::isalnum(data[offset]) || data[offset] == '_')) {
			str += data[offset++];
		}
		return Token{ Token::Kind::Identifier, str };
	}
}

std::vector<Token> tokenize(const std::string& data) {
	std::vector<Token> tokens;
	std::size_t offset = 0;
	while (offset < data.size()) {
		if (auto token = next_token(data, offset)) {
			tokens.push_back(std::move(*token));
		} else {
			break;
		}
	}
	return tokens;
}

vmf_group parse_vmf(std::string data) {
	auto tokens = tokenize(data);
	vmf_group root;
	std::stack<vmf_group*> groups;
	groups.push(&root);
	std::size_t ti = 0;
	while (ti < tokens.size()) {
		auto tok = tokens[ti];
		if (tok.kind == Token::Kind::Identifier) {
			auto& group = groups.top()->children.emplace(tok.value, vmf_group())->second;
			groups.push(&group);
			assert(tokens[ti + 1].kind == Token::Kind::LeftBrace);
			ti += 2;
			continue;
		} else if (tok.kind == Token::Kind::String) {
			groups.top()->kvs[tok.value] = tokens[ti + 1].value;
			ti += 2;
		} else if (tok.kind == Token::Kind::RightBrace) {
			groups.pop();
			++ti;
		} else {
			assert(false);
		}
	}
	return root;
}

vmf_group load_vmf(std::string path) {
	std::ifstream file(path.data());
	assert(file.good());
	std::string data, line;
	while (std::getline(file, line)) {
		data += line + "\n";
	}
	return parse_vmf(std::move(data));
}

// math

const double delta = 0.000001;
struct vec3 {
	double x, y, z;

	bool is_zero() const {
		return (std::abs(x) < delta && std::abs(y) < delta && std::abs(z) < delta);
	}

	double length() const {
		return std::sqrt(x * x + y * y + z * z);
	}
};

bool operator==(const vec3& a, const vec3& b) {
	return std::abs(a.x - b.x) < delta && std::abs(a.y - b.y) < delta && std::abs(a.z - b.z) < delta;
}

vec3 cross(const vec3& a, const vec3& b) {
	return vec3{
		a.y * b.z - a.z * b.y,
		a.z * b.x - a.x * b.z,
		a.x * b.y - a.y * b.x
	};
}

double dot(const vec3& a, const vec3& b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

vec3 operator-(const vec3& a, const vec3& b) {
	return vec3{ a.x - b.x, a.y - b.y, a.z - b.z };
}

bool parallel(const vec3& a, const vec3& b) {
	auto c = cross(a, b);
	return c.is_zero();
}

vec3 operator-(const vec3& a) {
	return { -a.x, -a.y, -a.z };
}

vec3 operator*(double a, const vec3& b) {
	return { b.x * a, b.y * a, b.z * a };
}

vec3 operator/(const vec3& a, double b) {
	return { a.x / b, a.y / b, a.z / b };
}

struct plane {
	vec3 points[3];
	std::string material;

	vec3 surface_normal() const {
		vec3 tri_edge_a = points[1] - points[0];
		vec3 tri_edge_b = points[2] - points[0];
		return cross(tri_edge_a, tri_edge_b);
	}
};

struct brush {
	std::vector<plane> planes;
};

// interpret vmf

plane parse_plane(std::string str) {
	plane p;
	std::sscanf(str.c_str(), "(%lf %lf %lf) (%lf %lf %lf) (%lf %lf %lf)",
		&p.points[0].x, &p.points[0].y, &p.points[0].z,
		&p.points[1].x, &p.points[1].y, &p.points[1].z,
		&p.points[2].x, &p.points[2].y, &p.points[2].z
	);
	return p;
}

std::vector<brush> load_brushes(const vmf_group vmf) {
	std::vector<brush> brushes;
	auto& world = vmf.children.find("world")->second;
	auto wcit = world.children.begin();
	while (wcit != world.children.end()) {
		if (wcit->first == "solid") {
			brush b;
			auto sit = wcit->second.children.begin();
			while (sit != wcit->second.children.end()) {
				if (sit->first == "side") {
					auto& kvs = sit->second.kvs;
					plane p = parse_plane(kvs.find("plane")->second);
					p.material = kvs.find("material")->second;
					b.planes.push_back(std::move(p));
				}
				++sit;
			}
			brushes.push_back(std::move(b));
		}
		++wcit;
	}
	return brushes;
}

// convert brush plane csg to face vertices

struct face {
	std::vector<vec3> points;
	std::string material;

	std::vector<std::pair<double, double>> uvs;
};

std::optional<vec3> three_plane_intersect(const plane& a, const plane& b, const plane& c) {
	auto an = a.surface_normal(), bn = b.surface_normal(), cn = c.surface_normal();
	if (std::abs(dot(an, cross(bn, cn))) < delta) {
		return std::nullopt;
	}
	auto ad = dot(an, a.points[0]), bd = dot(bn, b.points[0]), cd = dot(cn, c.points[0]);
	vec3 p = -ad * cross(bn, cn) - bd * cross(cn, an) - cd * cross(an, bn);
	p = p / dot(an, cross(bn, cn));
	return -p;
}

std::vector<vec3> unique_vec3s(const std::vector<vec3>& vecs) {
	std::vector<vec3> result;
	for (auto& vec : vecs) {
		bool unique = true;
		for (auto& result_vec : result) {
			if (result_vec == vec) {
				unique = false;
				break;
			}
		}
		if (unique) {
			result.push_back(vec);
		}
	}
	return result;
}

double angle_between(const vec3& a, const vec3& b, const vec3& normal) {
	auto result = std::acos(dot(a, b) / (a.length() * b.length()));
	auto d = dot(normal, cross(a, b));
	if (d < 0) {
		result *= -1;
	}
	return result;
}

std::vector<face> brush_to_faces(const brush& b) {
	std::vector<face> faces;
	for (auto& p1 : b.planes) {
		// ignore surfaces with tool materials
		if (p1.material.substr(0, 6) == "TOOLS/") {
			continue;
		}
		std::vector<vec3> vertices;
		// find all points where this plane intersects with two other planes from this brush
		for (auto& p2 : b.planes) {
			for (auto& p3 : b.planes) {
				if (auto intersection = three_plane_intersect(p1, p2, p3)) {
					vertices.push_back(*intersection);
				}
			}
		}
		// remote duplicate points (A,B,C vs A,C,B point)
		vertices = unique_vec3s(vertices);
		// remove any point on the surface normal side of any plane
		vertices.erase(std::remove_if(vertices.begin(), vertices.end(), [&](const vec3& v) {
			for (auto p : b.planes) {
				auto d = v - p.points[0];
				auto dist = dot(d, p.surface_normal());
				if (dist < 0 && std::abs(dist) > delta) { return true; }
			}
			return false;
		}), vertices.end());

		// build basis for 2d plane space
		auto nb_z = p1.surface_normal() / p1.surface_normal().length();
		vec3 center{ 0,0,0 };
		for (auto vert : vertices) {
			center.x += vert.x;
			center.y += vert.y;
			center.z += vert.z;
		}
		center = center / vertices.size();
		vec3 u = vertices[0] - center;
		u = u - dot(u, nb_z) * nb_z;
		auto nb_x = u / u.length();
		auto nb_y = cross(nb_z, nb_x);
		// sort vertices by angle from center point in plane space for winding order
		std::sort(vertices.begin(), vertices.end(), [&](const vec3& a, const vec3& b) {
			std::pair<double, double> a_proj{ dot(a - center, nb_x), dot(a - center, nb_y) };
			std::pair<double, double> b_proj{ dot(b - center, nb_x), dot(b - center, nb_y) };
			return std::atan2(a_proj.second, a_proj.first) > std::atan2(b_proj.second, b_proj.first);
		});

		// calculate some lame uvs
		std::vector<std::pair<double, double>> uvs;
		double min_u = INFINITY, max_u = -1 * INFINITY;
		double min_v = min_u, max_v = max_u;
		for (auto vert : vertices) {
			std::pair<double, double> proj{ dot(vert - center, nb_x), dot(vert - center, nb_y) };
			if (proj.first < min_u) { min_u = proj.first; }
			if (proj.second < min_v) { min_v = proj.second; }
			if (proj.first > max_u) { max_u = proj.first; }
			if (proj.second > max_v) { max_v = proj.second; }
		}
		double u_offset = -min_u;
		double u_range = max_u - min_u;
		double v_offset = -min_v;
		double v_range = max_v - min_v;
		for (auto vert : vertices) {
			std::pair<double, double> proj{ dot(vert - center, nb_x), dot(vert - center, nb_y) };
			double u = (proj.first + u_offset) / u_range;
			double v = (proj.second + v_offset) / v_range;
			uvs.push_back({ u, v });
		}

		face f;
		f.material = p1.material;
		f.points = vertices;
		f.uvs = uvs;
		faces.push_back(f);
	}
	return faces;
}

// write faces, vertices, and vertex uvs to file in .obj format
void write_obj(std::vector<face>& faces, std::string filename) {
	std::ofstream out(filename);
	assert(out.good());
	for (auto f : faces) {
		for (auto p : f.points) {
			out << "v " << p.x << " " << p.y << " " << p.z << "\n";
		}
		for (auto uv : f.uvs) {
			out << "vt " << uv.first << " " << uv.second << "\n";
		}
		out << "f";
		for (int i = f.points.size(); i != 0; --i) {
			out << " " << "-" << i;
			out << "/-" << i;
		}
		out << "\n";
	}
}

int main(int argc, char** argv) {
	if (argc != 3) {
		std::cerr << "usage: vmf2obj in.vmf out.obj\n";
		return EXIT_FAILURE;
	}
	std::string in_path = argv[1];
	std::string out_path = argv[2];
	auto g = load_vmf(in_path);
	auto b = load_brushes(g);
	std::vector<face> faces;
	for (auto bb : b) {
		auto fb = brush_to_faces(bb);
		faces.insert(faces.end(), fb.begin(), fb.end());
	}
	write_obj(faces, out_path);
	return EXIT_SUCCESS;
}
