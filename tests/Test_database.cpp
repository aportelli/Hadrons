#include <Hadrons/Database.hpp>

using namespace Grid;
using namespace Hadrons;

template <typename T>
void testType(const T &x, const std::string expected)
{
    LOG(Message) << x << ": " << SqlEntry::sqlType(x);
    if (SqlEntry::sqlType(x) == expected)
    {
        std::cout << " good" << std::endl;
    }
    else
    {
        std::cout << " bad" << std::endl;
        exit(EXIT_FAILURE);
    }
}

GRID_SERIALIZABLE_ENUM(TestEnum, undef, red, 1, blue, 2, green, 3);

struct TestStruct: Serializable
{
    GRID_SERIALIZABLE_CLASS_MEMBERS(TestStruct,
                                    TestEnum, e,
                                    int, x,
                                    double, y,
                                    bool , b);
};

struct TestEntry: SqlEntry
{
    HADRONS_SQL_FIELDS(TestEntry,
                       int, a,
                       float, b,
                       std::string, msg,
                       std::vector<int>, vec,
                       TestStruct, st);
};

int main(int argc, char *argv[])
{
    Grid_init(&argc, &argv);

    // test SQL type detection
    double             x = 1.567;
    float              y = 2.498;
    unsigned int       z = 34;
    std::string        s = "hello";
    std::vector<float> v = {0.1, 0.2, 0.3};

    testType(x, "REAL");
    testType(y, "REAL");
    testType(z, "INTEGER");
    testType(s, "TEXT");
    testType(v, "TEXT");

    // test SQL schema serialization
    TestEntry  entry;
    TestStruct st;

    st.e      = TestEnum::red;
    st.x      = 4;
    st.y      = 5.67;
    st.b      = true;
    entry.a   = 1;
    entry.b   = 2.45;
    entry.msg = "hello";
    entry.vec = {1, 2, 4};
    entry.st  = st;
    LOG(Message) << TestEntry::sqlSchema() << std::endl;
    LOG(Message) << entry.sqlInsert() << std::endl;

    Grid_finalize();
    
    return EXIT_SUCCESS;
}