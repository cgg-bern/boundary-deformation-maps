//
// Created by Valentin Nigolian on 16.10.2025.
//

#include "QualityVector.hh"

QualityVector::QualityVector(){}

QualityVector::QualityVector(std::initializer_list<double> list){
    for(auto val: list){
        insert_sort(val);
    }
}

int QualityVector::compare(const QualityVector& other,
                           double epsilon){

    auto& vec_a = *this;
    auto& vec_b = other;

    int smallest_vec(-2);
    int min_size(vec_a.size());
    if(vec_a.size() == vec_b.size()){
        smallest_vec = 0;
    }else if(vec_a.size() < vec_b.size()){
        min_size = vec_a.size();
        smallest_vec = 1;
    }else{
        min_size = vec_b.size();
        smallest_vec = -1;
    }
    //std::cout<<" -      epsilon: "<<epsilon<<std::endl;
    //std::cout<<" - smallest vec: "<<smallest_vec<<std::endl;

    for(int i(0); i < min_size; i++){
        double diff = vec_a[i] - vec_b[i];
        double rel_diff = diff / std::min(vec_a[i], vec_b[i]);
        //std::cout<<" -----------------"<<std::endl;
        //std::cout<<" diff (rel / abs) "<<i<<": "<<rel_diff<<" / "<<diff<<std::endl;
        if(std::abs(rel_diff) > epsilon){
            //std::cout<<" > eps"<<std::endl;
            if(diff < 0){
                return -1;
            }else{
                return 1;
            }
        }
        //std::cout<<" < eps"<<std::endl;
    }

    return smallest_vec;
}


void QualityVector::insert_sort(double value) {
    auto it = std::lower_bound(data_.begin(), data_.end(), value);
    data_.insert(it, value);
}

void QualityVector::merge(const QualityVector& other) {
    for (double value : other) {
        this->insert_sort(value);
    }
}


int QualityVector::size() const{
    return data_.size();
}

std::vector<double>::const_iterator QualityVector::cbegin() const {
    return data_.cbegin();
}

std::vector<double>::const_iterator QualityVector::cend() const {
    return data_.cend();
}

std::vector<double>::const_iterator QualityVector::begin() const {
    return data_.cbegin();
}

std::vector<double>::const_iterator QualityVector::end() const {
    return data_.cend();
}


const double& QualityVector::operator[](int i) const{
    return data_[i];
}

const double& QualityVector::back() const {
    return data_.back();
}


const double& QualityVector::front() const {
    return data_.front();
}

std::ostream& operator<<(std::ostream& os, const QualityVector& vec) {
    os << "[";
    for(int i(0); i<vec.size(); i++){
        os << vec[i];
        if(i != (vec.size()-1)){
            os<< ", ";
        }
    }
    for(auto val: vec){
    }
    os << "]";
    return os;
}