//
// Created by Valentin Nigolian on 16.10.2025.
//

#pragma once

#include "includes.hh"

class QualityVector {

public:
    QualityVector();

    QualityVector(std::initializer_list<double> list);

    /* Lexicographic comparison of two vectors sorted in ascending order.
 * @arg epsilon can be used to define how close two components can be before comparing the next one
 *
 * -1 -> this  < other,
 *  0 -> this == other, i.e. |a[i] - b[i]| < epsilon for all components
 *  1 -> this  > other
 *
 *  by convention, if one vector is LARGER than the other, and it's equal to the longer one, it's considered smaller
 *  e.g. {0,1,2} is smaller than {0,1,2,3}
     *  i.e., a vector of the same quality but with fewer elements is preferred
 */
    int compare(const QualityVector& other,
                double epsilon = DBL_EPSILON);

    void insert_sort(double value);

    void merge(const QualityVector& other);


    std::vector<double>::const_iterator cbegin() const;
    std::vector<double>::const_iterator cend() const;
    std::vector<double>::const_iterator begin() const;
    std::vector<double>::const_iterator end() const;

    int size() const;

    const double& operator[](int i) const;
    const double& back() const;
    const double& front() const;

private:
    std::vector<double> data_;
};


std::ostream& operator<<(std::ostream& os, const QualityVector& vec);