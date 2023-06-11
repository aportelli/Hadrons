/*
 * Test_database.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution 
 * directory.
 */

/*  END LEGAL */

#include <Hadrons/Database.hpp>

using namespace Grid;
using namespace Hadrons;

GRID_SERIALIZABLE_ENUM(TestEnum, undef, red, 1, blue, 2, green, 3);

struct TestStruct: Serializable
{
    GRID_SERIALIZABLE_CLASS_MEMBERS(TestStruct,
                                    TestEnum, e,
                                    int, x,
                                    double, y,
                                    bool , b);
};

struct IndexEntry: public SqlEntry
{
    HADRONS_SQL_FIELDS(unsigned int, traj);
};

struct TestEntry: public SqlEntry
{
    HADRONS_SQL_FIELDS(SqlNotNull<int>, a,
                       SqlNotNull<float>, b,
                       TestEnum, e,
                       std::string, msg,
                       SqlNotNull<std::vector<int>>, vec,
                       TestStruct, st);
};

#define TEST_TYPE(type, expected)\
LOG(Message) << #type << ": " << SqlEntry::sqlType<type>();\
if (SqlEntry::sqlType<type>() == expected)\
{\
    std::cout << " good" << std::endl;\
}\
else\
{\
    std::cout << " bad" << std::endl;\
    exit(EXIT_FAILURE);\
}\

int main(int argc, char *argv[])
{
    Grid_init(&argc, &argv);

    // test SQL type generation ////////////////////////////////////////////////
    TEST_TYPE(double, "REAL");
    TEST_TYPE(float, "REAL");
    TEST_TYPE(int, "INTEGER");
    TEST_TYPE(std::string, "TEXT");
    TEST_TYPE(std::vector<float>, "TEXT");
    TEST_TYPE(SqlNotNull<double>, "REAL NOT NULL");
    TEST_TYPE(SqlUnique<SqlNotNull<std::vector<double>>>, "TEXT NOT NULL UNIQUE");
    
    // test SQL schema serialization ///////////////////////////////////////////
    TestEntry  entry;
    TestStruct st;

    CppType<SqlNotNull<double>>::type x;

    st.e      = TestEnum::red;
    st.x      = 4;
    st.y      = 5.67;
    st.b      = true;
    entry.a   = 1;
    entry.b   = 2.45;
    entry.e   = TestEnum::blue;
    entry.msg = "hello";
    entry.vec = {1, 2, 4};
    entry.st  = st;
    LOG(Message) << TestEntry::sqlSchema() << std::endl;
    LOG(Message) << entry.sqlInsert() << std::endl;
    entry.nullify.e = true;
    LOG(Message) << entry.sqlInsert() << std::endl;
    entry.nullify.e = false;

    // test merging SqlEntry classes ///////////////////////////////////////////
    typedef MergedSqlEntry<IndexEntry, TestEntry> TableEntry;

    IndexEntry ind; ind.traj = 2000;
    TableEntry me = mergeSqlEntries(ind, entry);

    LOG(Message) << "schema               : " << me.sqlSchema() << std::endl;
    LOG(Message) << "insert               : " << me.sqlInsert() << std::endl;
    LOG(Message) << "cols                 : " << me.cols() << std::endl;
    LOG(Message) << "first entry (static) : " << me.getEntry<0>().sqlInsert() << std::endl;
    LOG(Message) << "first entry (dynamic): " << me.getEntry(0)->sqlInsert() << std::endl;

    // test Database class low-level SQL operations ////////////////////////////
    Database db("test.db");

    db.execute(
        "CREATE TABLE IF NOT EXISTS 'test' ("
        "   'name'	TEXT NOT NULL UNIQUE,   "
        "   'value' INTEGER NOT NULL        "
        ");");

    for (unsigned int i = 0; i < 10; ++i)
    {
        std::string is = std::to_string(i);

        db.execute("INSERT INTO 'test' VALUES ('i" + is + "', " + is + ");");
    }

    QueryResult table = db.execute("SELECT * FROM 'test'");
    std::string head;

    for (unsigned int j = 0; j < table.cols(); ++j)
    {
        head += table.colName(j) + "|";
    }
    head.pop_back();
    LOG(Message) << "Columns " << head << std::endl;
    for (unsigned int i = 0; i < table.rows(); ++i)
    {
        std::string msg;

        for (auto &col: table[i])
        {
            msg += col + "|";
        }
        msg.pop_back();
        LOG(Message) << msg << std::endl;
    }

    // test Database class high-level operations ///////////////////////////////
    db.createTable<TestEntry>("test2");
    for (unsigned int t = 1000; t < 2000; t += 20)
    {
        entry.msg  = "result_" + std::to_string(t);
        entry.st.x = t*2;
        if (t > 1200)
        {
            entry.nullify.e = true;
        }
        db.insert("test2", entry);
    }
    auto table2 = db.getTable<TestEntry>("test2");
    for (auto &e: table2)
    {
        LOG(Message) << e << std::endl;
    }
    auto col = db.getTableColumn<TestStruct>("test2", "st");
    for (auto &s: col)
    {
        LOG(Message) << s << std::endl;
    }

    LOG(Message) << "Table 'test' exists: " << db.tableExists("test") << std::endl;
    LOG(Message) << "Table 'foo' exists : " << db.tableExists("foo")  << std::endl;

    db.createKeyValueTable("kvtest");
    db.insertValue("kvtest", "someKey", st);
    auto buf = db.getValue<TestStruct>("kvtest", "someKey"); 
    assert(buf == st);

    Grid_finalize();
    
    return EXIT_SUCCESS;
}
