#include <string>

using namespace std;

class Molecule{
    public:
    int natom;
    int charge;
    int *zval;
    string regist;
    double **geom;
    string point_group;

    void print_geom();
    void rotate(double phi);
    void translate(double x,double y, double z);
    double bond(int atom1, int atom2);
    double center_of_mass();

    Molecule();
    ~Molecule();
};
