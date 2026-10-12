#pragma once

#include <exception>

template<typename _entity_handle, typename _value_type>
struct QueueNode{
    QueueNode(_entity_handle handle, _value_type value) : handle_(handle), value_(value), next_(nullptr){}

    _entity_handle handle_;
    _value_type value_;
    QueueNode* next_;
};


//NOTE: the lower the priority value, the higher the node is in the queue.
template<typename _entity_handle, typename _value_type>
class EntityPriorityQueue
{
public:
    using Node = QueueNode<_entity_handle, _value_type>;


    ~EntityPriorityQueue(){
        while(front_ != nullptr){
            auto current_front = front_;
            front_ = current_front->next_;
            delete current_front;
        }
    }


    void push(_entity_handle handle,
              _value_type value){

        Node* new_node = new Node(handle, value);
        
        if(!front_){
            front_ = new_node;
        }else{
            if(front_->value_ > value){
                new_node->next_ = front_;
                front_ = new_node;
            }else{
                push_after(new_node, front_);
            }
        }
    }
    
    void update_priority(_entity_handle handle, _value_type new_value){
        
        //first, remove the node
        if(front_->handle_ == handle){
            front_ = front_->next_;
        }else{
            Node* current = front_;
            Node* next = front_->next_;
            
            while(next != nullptr){
                if(next->handle_ == handle){
                    auto to_delete = next;
                    current->next_ = next->next_;
                    delete next;
                    break;
                }
                current = next;
                next = next->next_;
            }
        }

        //then re-insert updated node
        push(handle, new_value);
        
    }


    void pop(){
        if(!front_){
            return;
        }
        Node* to_remove = front_;
        front_ = front_->next_;
        delete to_remove;
    }


    //returns the node with the LOWEST priority
    _entity_handle front(){
        if(!front_){
            throw std::runtime_error("called 'front()' on empty queue");
        }else{
            return front_->handle_;
        }
    }

    _value_type front_priority(){
        if(!front_){
            throw std::runtime_error("called 'front_priority()' on empty queue");
        }else{
            return front_->value_;
        }
    }


    //TODO: make this O(1)
    int size() const{
        if(!front_){
            return 0;
        }else{
            int i(1);
            Node* next = front_->next_;
            while(next != nullptr){
                next = next->next_;
                i++;
            }
            return i;
        }
    }


    bool empty() const{
        return !front_;
    }

    void clear(){
        while(!empty()){
            pop();
        }
    }

#if 0
    void print(){
        if(front_){
            std::cout<<" LIST: "<<front_->value_;
            Node* next = front_->next_;

            while(next != nullptr){
                std::cout<<" -> "<<next->value_;
                next = next->next_;
            }
            std::cout<<std::endl;

        }else{
            std::cout<<" LIST : NIL"<<std::endl;
        }
    }
#endif



private:

    Node* front_ = nullptr;

    //assumes that after_node->value_ < to_insert->value_
    void push_after(Node* to_insert, Node* after_node){

        if(!after_node){//for safety
            return;
        }else if(!after_node->next_){
            after_node->next_ = to_insert;
        }else if(to_insert->value_ < after_node->next_->value_){
            to_insert->next_ = after_node->next_;
            after_node->next_ = to_insert;
        }else{
            push_after(to_insert, after_node->next_);
        }
    }

};

