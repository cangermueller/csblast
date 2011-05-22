#include <gtest/gtest.h>

#include "cs.h"
#include "blosum_matrix.h"
#include "profile-inl.h"
#include "count_profile-inl.h"

namespace cs {


class ProfileTest : public testing::Test {

  protected:
    ProfileTest() : ran(Ran(0)) {}
  
    Profile<AA> GetRndProfile(size_t len = kLen) {
      Profile<AA> p(len);
      for (size_t i = 0; i < len; ++i) {
        for (size_t a = 0; a < AA::kSize; ++a) 
          p[i][a] = ran.doub();
      }
      return p;
    }

    bool IsEqual(const Profile<AA>& p, const Profile<AA>& q) {
      if (p.length() != q.length()) return false;
      for (size_t i = 0; i < p.length(); ++i) {
        for (size_t a = 0; a < AA::kSize; ++a) {
          if (p[i][a] != q[i][a]) return false;
        }
      }
      return true;
    }


 protected:
    Ran ran;

    static const size_t kLen   = 100;
    static const double kDelta = 1e-5;
};

TEST_F(ProfileTest, InitValue) {
  const double value = ran.doub();
  Profile<AA> p(kLen, value);
  for (size_t i = 0; i < kLen; ++i) {
    for (size_t a = 0; a < AA::kSize; ++a)
      ASSERT_EQ(value, p[i][a]);
  }
}

TEST_F(ProfileTest, InitArray) {
  size_t l = kLen * AA::kSizeAny;
  double v[l];
  for (size_t i = 0; i < l; ++i)
    v[i] = ran.doub();
  Profile<AA> p(kLen, v);
  for (size_t i = 0; i < l; ++i) {
    ASSERT_EQ(p[0][i], v[i]);
    ASSERT_EQ(p[i / AA::kSizeAny][i % AA::kSizeAny], v[i]);
  }
}

TEST_F(ProfileTest, InitProfile) {
  Profile<AA> p = GetRndProfile();
  Profile<AA> q(p);
  ASSERT_TRUE(IsEqual(p, q));
}

TEST_F(ProfileTest, AssignProfile) {
  Profile<AA> p = GetRndProfile();
  Profile<AA> q;
  q = p;
  ASSERT_TRUE(IsEqual(p, q));
}

TEST_F(ProfileTest, InitCountProfile) {
  CountProfile<AA> cp(kLen);
  for (size_t i = 0; i < kLen; ++i) {
    for (size_t a = 0; a < AA::kSize; ++a)
      cp.counts[i][a] = ran(100);
    cp.neff[i] = ran(100) + 1;
  }
  Profile<AA> p(cp);
  ASSERT_EQ(cp.length(), p.length());
  for (size_t i = 0; i < cp.length(); ++i) {
    for (size_t a = 0; a < AA::kSize; ++a)
      ASSERT_EQ(cp.counts[i][a] / cp.neff[i], p[i][a]);
  }
}

TEST_F(ProfileTest, Insert) {
  Profile<AA> dest(kLen);
  for (size_t r = 0; r < 10; ++r) {
    Profile<AA> src = GetRndProfile();
    size_t idx = ran(kLen);
    dest.Insert(idx, src);
    for (size_t i = idx; i < kLen; ++i) {
      for (size_t a = 0; a < AA::kSize; ++a)
        ASSERT_EQ(src[i - idx][a], dest[i][a]);
    }
  }
}

TEST_F(ProfileTest, Neff) {
  Profile<AA> p(kLen, 0.0);
  for (size_t i = 0; i < kLen; ++i)
    p[i][ran(AA::kSize)] = 1.0;
  EXPECT_NEAR(1.0, Neff(p), kDelta);
  Assign(p, 0.0);
  for (size_t i = 0; i < kLen; ++i) {
    size_t a = ran(AA::kSize);
    p[i][a] = 0.5;
    p[i][(a + 1) % AA::kSize] = 0.5;
  }
  EXPECT_NEAR(2.0, Neff(p), kDelta);
}





};  // namespace cs
