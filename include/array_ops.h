#include <boost/python.hpp>
void nmat_detvecs_apply(const bp::object & ft, const bp::object & bins, const bp::object & iD, const bp::object & iV, float s, float nsamp);
int process_cuts(const bp::object & range_matrix, const std::string & operation, const std::string & model, const bp::dict & params, const bp::object & tod, const bp::object & vals);
void test_cuts(const bp::object & range_matrix);
