#include "main_test.h"

TEST(CHECK, boolator_a) {
    ASSERT_EQ(true, boolator_a());
}
TEST(CHECK, boolator_b) {
    ASSERT_EQ(true, boolator_b());
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}

