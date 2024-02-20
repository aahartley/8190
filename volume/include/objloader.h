#ifndef OBJ_LOADER_H
#define OBJ_LOADER_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include "Vector.h"
#include "Models.h"

namespace lux
{

struct Triangle
{
  Triangle(){}
  Triangle(Vector vv1, Vector vv2, Vector vv3) : v1(vv1), v2(vv2), v3(vv3)
  {
    e1 = v2 - v1;
    e2 = v3 - v1;
    e3 = v3 - v2;
  }

  Vector v1, v2, v3, e1, e2, e3;
};
class ObjLoader {
public:
    static std::vector<Triangle> loadObj(const std::string& filename) {
        std::vector<Vector> vertices;
        std::vector<Triangle> triangles;

        std::ifstream file(filename);
        std::string line;

        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string token;
            iss >> token;

            if (token == "v") {
                Vector vertex;
                float x,y,z;
                iss >> x >> y >> z;
                vertex.set(x,y,z);
                vertices.push_back(vertex);
            } else if (token == "f") {
                int v1Index, v2Index, v3Index;
                iss >> v1Index  >> v2Index  >> v3Index;
                Vector v1 = vertices[v1Index - 1];
                Vector v2 = vertices[v2Index - 1];
                Vector v3 = vertices[v3Index - 1];
                triangles.push_back(Triangle(v1,v2,v3));
                //std::cout << v1Index << ' ' << v2Index << ' ' << v3Index << ' ' << v1.X() << ' ' << v1.Y() << ' ' << v1.Z() << '\n';
            }
        }

        file.close();
        return triangles;
    }
};
    
} // namespace lux
#endif // OBJ_LOADER_H