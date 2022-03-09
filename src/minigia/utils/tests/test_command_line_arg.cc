#include <catch2/catch_test_macros.hpp>

#include <catch2/catch_approx.hpp>
using Catch::Approx;

#include "../command_line_arg.hpp"

TEST_CASE("test command line argument parse", "[Command_line_arg]")
{

    Command_line_arg arg("foo");

    SECTION("get_lhs")
    {
        Command_line_arg arg("foo=bar");
        REQUIRE(arg.get_lhs() == "foo");
    }

    SECTION("get_rhs")
    {
        Command_line_arg arg("foo=bar");
        REQUIRE(arg.get_rhs() == "bar");
    }

    SECTION("extract_value_int")
    {
        Command_line_arg arg("foo=7");
        REQUIRE(arg.extract_value<int>() == 7);
    }

    SECTION("extract_value_double")
    {
        Command_line_arg arg("foo=2.718");

        REQUIRE(arg.extract_value<double>() == Approx(2.718).epsilon(1e-13));
    }

    SECTION("extract_value_bool")
    {
        Command_line_arg arg("foo=1");
        REQUIRE(arg.extract_value<bool>() == true);

        Command_line_arg arg2("foo=0");
        REQUIRE(arg2.extract_value<bool>() == false);
    }

    SECTION("extract_value_string")
    {
        Command_line_arg arg("foo=bar");
        std::string bar("bar");
        REQUIRE(arg.extract_value<std::string>() == bar);
    }
}
