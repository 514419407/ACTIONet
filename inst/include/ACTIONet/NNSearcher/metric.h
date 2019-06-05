#ifndef METRIC_H
#define METRIC_H

#include <algorithm>
#include <climits>

#include <nn_aux.h>

// This header file defines templated function objects (functors) for distance
// calculations. They work on any reasonable container class by using forward
// iterators For good performance, no bound checking is done, so it's the
// programmer's resposibility to have everything properly allocated.
//
// To support PARTIAL DISTANCE CALCULATION during the search phase (not the
// preprocessing phase!), operator() is overloaded. operator() with three
// arguments is the standard distance calculation and returns the exact
// distance. operator() with four arguments is the tresholded distance
// calculation with terminates as soon as the partial distance exceeds the given
// threshold value (but otherwise returns the exact distance). Carefull
// implementation can speed up the search considerably. When threshold is
// exceeded, DBL_MAX is returned to prevent the search functions to consider
// this value as valid distance.


class jensen_distance {
public:
  jensen_distance(){};
  template <class ForwardIterator1, class ForwardIterator2>
  double operator()(ForwardIterator1 first1, ForwardIterator1 last1,
                  ForwardIterator2 first2) const {
	
	int idx = 0;
	
	double mid = 0.0, p_div = 0.0, q_div = 0.0;	
    for (++first1, ++first2; first1 != last1; ++first1, ++first2) {
		
		double p = *first1;
		double q = *first2;
		
		double mid = (p + q) / 2;		
		
        if(mid != 0 && p != 0) {
			p_div += p*std::log2(p / mid);
        }
        
        if(mid != 0 && q != 0) {
			q_div += q*std::log2(q / mid);        
		}		
    }
    
    double JS = (p_div + q_div) / 2.0;
	
    return sqrt(JS); // To make it a "metric"
  }

  // support partial search : if partial distance exceeds thresh, stop computing
  // of distance
  template <class ForwardIterator1, class ForwardIterator2>
  double operator()(ForwardIterator1 first1, ForwardIterator1 last1,
                  ForwardIterator2 first2, const double thresh) const {
	
	double t = 2.0 * (thresh * thresh);
	
	double mid = 0.0, sum = 0.0;
    for (++first1, ++first2; first1 != last1; ++first1, ++first2) {
		double p = *first1;
		double q = *first2;
		
		double mid = (p + q) / 2;		
		
        if(mid != 0 && p != 0) {
			sum += p*std::log2(p / mid);
        }
        
        if(mid != 0 && q != 0) {
			sum += q*std::log2(q / mid);        
		}		
		
		if (sum > t)
		  return DBL_MAX;
	}
    
    return( sqrt(sum / 2.0) );
  }
};


class euclidian_distance {
public:
  euclidian_distance(){};
  template <class ForwardIterator1, class ForwardIterator2>
  double operator()(ForwardIterator1 first1, ForwardIterator1 last1,
                  ForwardIterator2 first2) const {
    const double y = (*first1 - *first2);
    double dist = y * y;
    for (++first1, ++first2; first1 != last1; ++first1, ++first2) {
      const double x = (*first1 - *first2);
      dist += x * x;
    }
    return sqrt(dist);
  }

  // support partial search : if partial distance exceeds thresh, stop computing
  // of distance
  template <class ForwardIterator1, class ForwardIterator2>
  double operator()(ForwardIterator1 first1, ForwardIterator1 last1,
                  ForwardIterator2 first2, const double thresh) const {
    const double t = thresh * thresh;
    const double y = (*first1 - *first2);
    double dist = y * y;

    if (dist > t)
      return DBL_MAX;

    for (++first1, ++first2; first1 != last1; ++first1, ++first2) {
      const double x = (*first1 - *first2);
      dist += x * x;
      if (dist > t)
        return DBL_MAX;
    }
    return sqrt(dist);
  }
};

class maximum_distance {
public:
  maximum_distance(){};
  template <class ForwardIterator1, class ForwardIterator2>
  double operator()(ForwardIterator1 first1, ForwardIterator1 last1,
                    ForwardIterator2 first2) const {
    double dist = fabs(*first1 - *first2);
    for (++first1, ++first2; first1 != last1; ++first1, ++first2) {
      const double x = fabs(*first1 - *first2);
      if (x > dist)
        dist = x;
    }
    return dist;
  }
  // support partial search
  template <class ForwardIterator1, class ForwardIterator2>
  double operator()(ForwardIterator1 first1, ForwardIterator1 last1,
                    ForwardIterator2 first2, const double thresh) const {
    double dist = fabs(*first1 - *first2);
    if (dist > thresh) {
      return DBL_MAX;
    }
    for (++first1, ++first2; first1 != last1; ++first1, ++first2) {
      const double x = fabs(*first1 - *first2);
      if (x > dist) {
        if (x > thresh) {
          return DBL_MAX;
        }
        dist = x;
      }
    }
    return dist;
  }
};

class manhattan_distance {
public:
  manhattan_distance(){};
  template <class ForwardIterator1, class ForwardIterator2>
  double operator()(ForwardIterator1 first1, ForwardIterator1 last1,
                    ForwardIterator2 first2) const {
    double dist = fabs(*first1 - *first2);
    for (++first1, ++first2; first1 != last1; ++first1, ++first2) {
      const double x = fabs(*first1 - *first2);
      dist += x;
    }
    return dist;
  }
  // support partial search
  template <class ForwardIterator1, class ForwardIterator2>
  double operator()(ForwardIterator1 first1, ForwardIterator1 last1,
                    ForwardIterator2 first2, const double thresh) const {
    double dist = fabs(*first1 - *first2);
    if (dist > thresh) {
      return DBL_MAX;
    }
    for (++first1, ++first2; first1 != last1; ++first1, ++first2) {
      const double x = fabs(*first1 - *first2);
      dist += x;
      if (dist > thresh) {
        return DBL_MAX;
      }
    }
    return dist;
  }
};

class hamming_distance {
public:
  hamming_distance(){};
  template <class ForwardIterator1, class ForwardIterator2>
  double operator()(ForwardIterator1 first1, ForwardIterator1 last1,
                    ForwardIterator2 first2) const {
    int dist = 0;
    for (; first1 != last1; ++first1, ++first2) {
      if (*first1 != *first2) {
        dist++;
      }
    }
    return dist;
  }
  // support partial search
  template <class ForwardIterator1, class ForwardIterator2>
  double operator()(ForwardIterator1 first1, ForwardIterator1 last1,
                    ForwardIterator2 first2, const double thresh) const {
    int dist = 0;
    for (; first1 != last1; ++first1, ++first2) {
      if (*first1 != *first2) {
        dist++;
        if (dist > thresh) {
          return DBL_MAX;
        }
      }
    }
    return dist;
  }
};

class euclidian_distance_unrolled {
public:
  euclidian_distance_unrolled(){};
  template <class ForwardIterator1, class ForwardIterator2>
  double operator()(ForwardIterator1 first1, const ForwardIterator1 last1,
                  ForwardIterator2 first2) const {
    double dist = 0, diff = 0;

    size_t n4 = (last1-first1)/4;
    const ForwardIterator1 partial1 = first1 + n4 * 4;

    while (first1 != partial1) {
      diff = *first1 - *first2; dist += diff * diff;
      ++first1; ++first2;
      diff = *first1 - *first2; dist += diff * diff;
      ++first1; ++first2;
      diff = *first1 - *first2; dist += diff * diff;
      ++first1; ++first2;
      diff = *first1 - *first2; dist += diff * diff;
      ++first1; ++first2;
    }

    for (; first1 != last1; ++first1, ++first2) {
      diff = (*first1 - *first2);
      dist += diff * diff;
    }
    return sqrt(dist);
  }

  // Supports partial search : if partial distance exceeds thresh, stop computing
  // of distance
  template <class ForwardIterator1, class ForwardIterator2>
  double operator()(ForwardIterator1 first1, const ForwardIterator1 last1,
                  ForwardIterator2 first2, const double thresh) const {
    const double t = thresh * thresh;
    double dist = 0, diff = 0;

    size_t n4 = (last1-first1)/4;
    const ForwardIterator1 partial1 = first1 + n4 * 4;

    while (first1 != partial1) {
      diff = *first1 - *first2; dist += diff * diff;
      ++first1; ++first2;
      diff = *first1 - *first2; dist += diff * diff;
      ++first1; ++first2;
      diff = *first1 - *first2; dist += diff * diff;
      ++first1; ++first2;
      diff = *first1 - *first2; dist += diff * diff;
      ++first1; ++first2;
      if (dist > t)
        return DBL_MAX;
    }
    for (; first1 != last1; ++first1, ++first2) {
      diff = (*first1 - *first2);
      dist += diff * diff;
    }
    return sqrt(dist);
  }
};

#endif
