#include <pybindings.h>

#include <iostream>
#include <boost/python.hpp>

#include <container_pybindings.h>
#include <hkagg.h>
#include <exceptions.h>

// Templates for vector operations.  g3vectype might be, for example,
// G3VectorDouble.  One could go all the way and SFINAE this... but
// perhaps not everything should require a black belt.

template <typename g3vectype>
inline
int g3_vect_size(const G3FrameObjectPtr vp)
{
    auto v = boost::dynamic_pointer_cast<const g3vectype>(vp);
    return v == nullptr ? -1 : v->size();
}

template <typename g3vectype>
inline
void g3_concat(g3vectype &output, const g3vectype &src1, const g3vectype &src2)
{
    auto dest = output.begin();
    for (auto p1: src1)
        *(dest++) = p1;
    for (auto p2: src2)
        *(dest++) = p2;
}

template <typename g3vectype>
inline
G3FrameObjectPtr test_and_concat(const G3FrameObjectPtr src1, const G3FrameObjectPtr src2)
{
    auto v1 = boost::dynamic_pointer_cast<const g3vectype>(src1);
    auto v2 = boost::dynamic_pointer_cast<const g3vectype>(src2);
    if (v1 == nullptr || v2 == nullptr)
        return nullptr;
    boost::shared_ptr<g3vectype> outputp(new g3vectype());
    outputp->resize(v1->size() + v2->size());
    g3_concat(*outputp, *v1, *v2);
    return outputp;
}




/* IrregBlock */

std::string IrregBlock::Description() const
{
	std::ostringstream s;
        s << "<co-sampled vectors with " << times.size() << " samples>{";
	for (auto i = this->begin(); i != this->end(); ) {
            s << i->first;
            if (++i != this->end())
                s << ", ";
        }
        s << "}";
	return s.str();
}

std::string IrregBlock::Summary() const
{
    return Description();
}

template <class A> void IrregBlock::save(A &ar, unsigned v) const
{
	using namespace cereal;
	ar & make_nvp("parent", base_class<G3MapFrameObject>(this));
	ar & make_nvp("times", times);
}
template <class A> void IrregBlock::load(A &ar, unsigned v)
{
	using namespace cereal;
	ar & make_nvp("parent", base_class<G3MapFrameObject>(this));
	ar & make_nvp("times", times);
}

bool IrregBlock::Check()
{
    int n = times.size();
    
    // How to polymorph?  This is how.
    for (auto item = begin(); item != end(); ++item) {
        auto name = item->first;
        auto el = item->second;

        int check_type = 0;
        int check_len = -1;

        // Try to get a length...
        if ((check_len = g3_vect_size<G3VectorDouble>(el)) >= 0 ||
            (check_len = g3_vect_size<G3VectorInt   >(el)) >= 0 ||
            (check_len = g3_vect_size<G3VectorString>(el)) >= 0)
            check_type = 1;

        if (check_len >= 0 && check_len != n) {
            std::ostringstream s;
            s << "Vector not same length as .times: " << name << "\n";
            throw general_agreement_exception(s.str());
        }

        if (check_type) {
            std::ostringstream s;
            s << "Item type is not acceptable: " << name << "\n";
            throw general_type_exception(s.str());
        }
    }
    return true;
}

const
IrregBlock IrregBlock::Concatenate(const IrregBlock other)
{
    // Check that all keys in other are in this.
    for (auto item = other.begin(); item != other.end(); ++item) {
        if (find(item->first) == end()) {
            std::ostringstream s;
            s << "inconsistent keys: " << item->first << " on right only.";
            throw general_agreement_exception(s.str());
        }
    }
    
    int n_cat = times.size() + other.times.size();
    IrregBlock output;
    output.times.resize(n_cat);
    g3_concat(output.times, times, other.times);

    for (auto item = begin(); item != end(); ++item) {
        auto oitem = other.find(item->first);
        if (oitem == other.end()) {
            std::ostringstream s;
            s << "inconsistent keys: " << item->first << " on left only.";
            throw general_agreement_exception(s.str());
        }

        G3FrameObjectPtr catted;
        if (
            (catted = test_and_concat<G3VectorDouble>(
                item->second, oitem->second)) != nullptr ||
            (catted = test_and_concat<G3VectorInt>(
                item->second, oitem->second)) != nullptr ||
            (catted = test_and_concat<G3VectorString>(
                item->second, oitem->second)) != nullptr
            ) {
            output.insert(std::make_pair(item->first, catted));
        } else {
            std::ostringstream s;
            s << "unhandled type for key: " << item->first;
            throw general_agreement_exception(s.str());
        }
    }

    return output;
}


/* IrregBlockDouble */

std::string IrregBlockDouble::Description() const
{
	std::ostringstream s;
	s << "Double data (" << data.size() << " vectors) with timestamp.";
	return s.str();
}

std::string IrregBlockDouble::Summary() const
{
    return Description();
}

template <class A> void IrregBlockDouble::serialize(A &ar, unsigned v)
{
	using namespace cereal;
        // v is the version code!

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("prefix", prefix);
	ar & make_nvp("t", t);
	ar & make_nvp("data", data);
}


G3_SPLIT_SERIALIZABLE_CODE(IrregBlock);
G3_SERIALIZABLE_CODE(IrregBlockDouble);


namespace bp = boost::python;

PYBINDINGS("so3g")
{
    // This is based on register_g3map macro.
    bp::class_<IrregBlock, bp::bases<G3FrameObject,
      std::map<typename IrregBlock::key_type, typename IrregBlock::mapped_type> >,
           boost::shared_ptr<IrregBlock> >("IrregBlock", "docstring")
        .def(bp::init<const IrregBlock &>())
        .def(bp::std_map_indexing_suite<IrregBlock, true>())
        .def_pickle(g3frameobject_picklesuite<IrregBlock>())
        // Extensions for IrregBlock are here:
        .def_readwrite("times", &IrregBlock::times, "Timestamp vector.")
        .def("Check", &IrregBlock::Check, "Check for internal consistency.")
        .def("Concatenate", &IrregBlock::Concatenate, "Concatenate two compatible IrregBlock.")
        ;
    register_pointer_conversions<IrregBlock>();


    EXPORT_FRAMEOBJECT(IrregBlockDouble, init<>(),
    "Data block for irregularly sampled data.")
    .def_readwrite("prefix", &IrregBlockDouble::prefix,
    "Prefix for field names.")
    .def_readwrite("data", &IrregBlockDouble::data,
    "Map to HK data vectors.")
    .def_readwrite("t", &IrregBlockDouble::t,
    "Timestamp vector.")
    ;

    bp::enum_<HKFrameType>("HKFrameType",
                           "Identifier for generic HK streams.")
        .value("session", HKFrameType::session)
        .value("status",  HKFrameType::status)
        .value("data",    HKFrameType::data)
        ;

}
