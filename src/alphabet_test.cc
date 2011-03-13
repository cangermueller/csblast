#include <gtest/gtest.h>

#include "cs.h"
#include "aa.h"

namespace cs {

TEST(AlphabetTest, AminoAcidAlphabet) {
  const char chars[] = {
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X'
  };

  ASSERT_EQ(static_cast<size_t>(20), AA::kSize);
  ASSERT_EQ(static_cast<size_t>(21), AA::kSizeAny);
  for (size_t a = 0; a < AA::kSizeAny; ++a)
    EXPECT_EQ(chars[a], AA::kIntToChar[AA::kCharToInt[chars[a]]]);
  for (size_t a = 0; a < AA::kSizeAny; ++a)
    EXPECT_EQ(a, AA::kCharToInt[AA::kIntToChar[a]]);
}

TEST(AlphabetTest, AS62) {
  ASSERT_EQ(static_cast<size_t>(62), AS62::kSize);
  ASSERT_EQ(static_cast<size_t>(63), AS62::kSizeAny);
  for (size_t a = 0; a < AS62::kEndGap; ++a) {
    LOG(ERROR) << a << "\t" << AS62::kIntToChar[a] << "\t" << (int) AS62::kCharToInt[(int) AS62::kIntToChar[a]];
    EXPECT_EQ(a, static_cast<size_t>(AS62::kCharToInt[AS62::kIntToChar[a]]));
  }
  for (size_t a = 0; a < AS62::kSizeAny; ++a)
    EXPECT_EQ(AS62::kIntToChar[a], AS62::kIntToChar[AS62::kCharToInt[AS62::kIntToChar[a]]]);
}

TEST(AlphabetTest, AS90) {
  ASSERT_EQ(static_cast<size_t>(90), AS90::kSize);
  ASSERT_EQ(static_cast<size_t>(91), AS90::kSizeAny);
  for (size_t a = 0; a < AS90::kSize; ++a) {
    LOG(ERROR) << a << "\t" << (int) AS90::kIntToChar[a] << "\t" << AS90::kIntToChar[a] << "\t" << (int) AS90::kCharToInt[(int) AS90::kIntToChar[a]];
    EXPECT_EQ(a, static_cast<size_t>(AS90::kCharToInt[AS90::kIntToChar[a]]));
  }
  for (size_t a = 0; a < AS90::kSizeAny; ++a)
    EXPECT_EQ(AS90::kIntToChar[a], AS90::kIntToChar[AS90::kCharToInt[AS90::kIntToChar[a]]]);
}

}  // namespace cs
