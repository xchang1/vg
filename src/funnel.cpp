#include "funnel.hpp"

#include <cassert>
#include <cstring>

/**
 * \file funnel.hpp: implementation of the Funnel class
 */
 
namespace vg {
using namespace std;

void Funnel::start(const string& name) {
    assert(!name.empty());

    // (Re)start the funnel.
    funnel_name = name;
    
    // Clear out old data
    stage_name.clear();
    substage_name.clear();
    stages.clear();
}

void Funnel::stop() {
    // Stop any lingering stage (which takes care of substages)
    stage_stop();
}

void Funnel::stage(const string& name) {
    assert(!funnel_name.empty());
    assert(!name.empty());
    
    // Stop the previous stage if any.
    stage_stop();

    // Allocate new stage structures.
    stages.emplace_back();
    stages.back().name = name;
    
    // Save the name
    stage_name = name;
}

void Funnel::stage_stop() {
    if (!stage_name.empty()) {
        // A stage was running.
        
        // Stop any substage.
        substage_stop();
        // Stop any process/produce 
        processed_input();
        produced_output();
        
        // Say the stage is stopped 
        stage_name.clear();
    }
}

void Funnel::substage(const string& name) {
    assert(!funnel_name.empty());
    assert(!stage_name.empty());
    assert(!name.empty());

    // Stop previous substage if any.
    substage_stop();
    
    // Substages don't bound produce/process.
    
    // Save the name 
    substage_name = name;
}
    
void Funnel::substage_stop() {
    if (!substage_name.empty()) {
        // A substage was running.
        
        // Substages don't bound produce/process.
        
        // Say the stage is stopped 
        substage_name.clear();
    }
}

void Funnel::processing_input(size_t prev_stage_item) {
    // We can only take input from previous stages, in a stage
    assert(!stage_name.empty());
    assert(stages.size() > 1);
    assert(prev_stage_item != numeric_limits<size_t>::max());
    assert(stages[stages.size() - 2].items.size() > prev_stage_item);

    // Stop any previous input processing
    processed_input();
    
    // Start this one
    input_in_progress = prev_stage_item;
}

void Funnel::processed_input() {
    if (input_in_progress != numeric_limits<size_t>::max()) {
        // We were processing an input
        
        // Say we're done with the input.
        input_in_progress = numeric_limits<size_t>::max();
    }
}

void Funnel::producing_output(size_t item) {
    // We can only produce output in a stage
    assert(!stage_name.empty());
    assert(!stages.empty());
    assert(item != numeric_limits<size_t>::max());
    
    // Stop any previous input processing
    produced_output();
    
    // Start this one
    output_in_progress = item;
}

void Funnel::produced_output() {
    if (output_in_progress != numeric_limits<size_t>::max()) {
        // We were producing an output
        
        // Say we're done with the output.
        output_in_progress = numeric_limits<size_t>::max();
    }
}

void Funnel::introduce(size_t count) {
    // Create that many new items
    for (size_t i = 0; i < count; i++) {
        create_item();
    }
}

void Funnel::expand(size_t prev_stage_item, size_t count) {
    for (size_t i = 0; i < count; i++) {
        // Create the requested number of items
        project(prev_stage_item);
    }
}

void Funnel::project(size_t prev_stage_item) {
    // There must be a prev stage to project from
    assert(stages.size() > 1);
    auto& prev_stage = stages[stages.size() - 2];

    // Make one new item
    size_t index = create_item();

    // Record the ancestry
    get_item(index).prev_stage_items.push_back(prev_stage_item);

    if (prev_stage.items[prev_stage_item].correct) {
        // Tag the new item correct if it came from something correct
        tag_correct(index);
    }
}

void Funnel::project_group(size_t prev_stage_item, size_t group_size) {
    // Project the item
    project(prev_stage_item);
    // Save the group size
    get_item(latest()).group_size = group_size;
}

void Funnel::fail(const char* filter, size_t prev_stage_item, double statistic) {
    // There must be a prev stage to project from
    assert(stages.size() > 1);
    auto& prev_stage = stages[stages.size() - 2];

    // Record the item as having failed this filter
    prev_stage.items[prev_stage_item].failed_filter = filter;
    prev_stage.items[prev_stage_item].failed_statistic = statistic;
}

void Funnel::pass(const char* filter, size_t prev_stage_item, double statistic) {
    // There must be a prev stage to project from
    assert(stages.size() > 1);
    auto& prev_stage = stages[stages.size() - 2];

    // Record the item as having passed this filter
    prev_stage.items[prev_stage_item].passed_filters.emplace_back(filter);
    prev_stage.items[prev_stage_item].passed_statistics.emplace_back(statistic);
}

void Funnel::score(size_t item, double score) {
    get_item(item).score = score;
}

void Funnel::tag_correct(size_t item) {
    // Say the item is correct
    get_item(item).correct = true;
    // Say the stage has something correct.
    stages.back().has_correct = true;
}

string Funnel::last_correct_stage() const {
    // Just do a linear scan backward through stages
    for (auto it = stages.rbegin(); it != stages.rend(); ++it) {
        if (it->has_correct) {
            return it->name;
        }
    }
    return "none";
}

size_t Funnel::latest() const {
    assert(!stages.empty());
    assert(!stages.back().items.empty());
    return stages.back().items.size() - 1;
}

void Funnel::for_each_stage(const function<void(const string&, const vector<size_t>&)>& callback) const {
    for (auto& stage : stages) {
        // Make a vector of item sizes
        vector<size_t> item_sizes;
        item_sizes.reserve(stage.items.size());
        for (auto& item : stage.items) {
            item_sizes.push_back(item.group_size);
        }
        // Report the name and item count of each stage.
        callback(stage.name, item_sizes);
    }
}

void Funnel::for_each_filter(const function<void(const string&, const string&,
    const FilterPerformance&, const FilterPerformance&, const vector<double>&, const vector<double>&)>& callback) const {
    
    for (auto& stage : stages) {
        // Hold the names of all filters encountered
        vector<const char*> filter_names;
        // And the by-item and by-size performance stats.
        vector<pair<FilterPerformance, FilterPerformance>> filter_performances;
        // And the correct and not-known-correct filter statistic values
        vector<pair<vector<double>, vector<double>>> filter_statistics;
        
        for (auto& item : stage.items) {
            // For each item
            size_t filter_index;
            for (filter_index = 0; filter_index < item.passed_filters.size(); filter_index++) {
                // For each filter it passed
                if (filter_index >= filter_names.size()) {
                    // If it is new
                
                    // Remember its name in the list of filters
                    filter_names.push_back(item.passed_filters[filter_index]);
                    
                    // And give it an empty report
                    filter_performances.emplace_back();
                    filter_statistics.emplace_back();
                } else {
                    // Make sure the name is correct
                    // TODO: can we justy match on pointer value and not string value?
                    assert(strcmp(filter_names[filter_index], item.passed_filters[filter_index]) == 0);
                }
                
                // Record passing
                filter_performances[filter_index].first.passing++;
                filter_performances[filter_index].first.passing_correct += item.correct;
                
                filter_performances[filter_index].second.passing += item.group_size;
                filter_performances[filter_index].second.passing_correct += item.correct ? item.group_size : 0;
                
                if (item.correct) {
                    // Record this statistic value as belonging to a correct item
                    filter_statistics[filter_index].first.push_back(item.passed_statistics[filter_index]);
                } else {
                    // Record this statistic value as belonging to a not necessarily correct item
                    filter_statistics[filter_index].second.push_back(item.passed_statistics[filter_index]);
                }
            }
            
            if (item.failed_filter != nullptr) {
                // For the final, failed filter, if any
                
                if (filter_index >= filter_names.size()) {
                    // If it is new
                    
                    // Remember its name in the list of filters
                    filter_names.push_back(item.failed_filter);
                    
                    // And give it an empty report
                    filter_performances.emplace_back();
                    filter_statistics.emplace_back();
                } else {
                    // Make sure the name is correct
                    // TODO: can we justy match on pointer value and not string value?
                    assert(strcmp(filter_names[filter_index], item.failed_filter) == 0);
                }
                
                // Record failing
                filter_performances[filter_index].first.failing++;
                filter_performances[filter_index].first.failing_correct += item.correct;
                
                filter_performances[filter_index].second.failing += item.group_size;
                filter_performances[filter_index].second.failing_correct += item.correct ? item.group_size : 0;
                
                if (item.correct) {
                    // Record this statistic value as belonging to a correct item
                    filter_statistics[filter_index].first.push_back(item.failed_statistic);
                } else {
                    // Record this statistic value as belonging to a not necessarily correct item
                    filter_statistics[filter_index].second.push_back(item.failed_statistic);
                }
            }
        }
        
        // Now we have gone through the filters for this stage for every item.
        
        for (size_t i = 0; i < filter_names.size(); i++) {
            // For each filter
            
            // Report the results tabulated across items.
            callback(stage.name, filter_names[i],
                filter_performances[i].first, filter_performances[i].second,
                filter_statistics[i].first, filter_statistics[i].second);
        }
    }
}

void Funnel::to_dot(ostream& out) {
    out << "digraph graphname {" << endl;
    out << "rankdir=\"TB\";" << endl;

    for (size_t s = 0; s < stages.size(); s++) {
        // For each stage in order
        auto& stage = stages[s];

        // Compute a GraphViz ID part for the stage
        string stage_id = "s" + to_string(s);

        // Start a subgraph.
        // Prepend cluster so it draws as a box.
        out << "subgraph cluster_" << stage_id << " {" << endl;
        out << "label = \"" << stage.name << "\";" << endl;
        out << "graph[style=solid];" << endl;
        out << "rank=same;" << endl;

        for (size_t i = 0; i < stage.items.size(); i++) {
            // For each item in the stage
            auto& item = stage.items[i];

            // Compute a GraphViz ID
            string item_id = stage_id + "i" + to_string(i);

            // Emit a node
            out << item_id << "[label=\"" << i << "\" shape=circle tooltip=\"";
            if (item.group_size != 0) {
                out << "size " << item.group_size;
            }
            if (item.score != 0) {
                if (item.group_size != 0) {
                    out << ", ";
                }
                out << "score " << item.score;
            }
            out << "\"";
            if (item.correct) {
                // Make it green if it is correct
                out << " color=green";
            }
            out << "];" << endl;

            if (s > 0) {
                // There is a previous stage, so we can draw edges from it.
                for (auto& p : item.prev_stage_items) {
                    // Connect everything from the previous stage to it
                    auto& prev_item = stages[s - 1].items.at(p);

                    out << "s" << (s - 1) << "i" << p << " -> " << item_id << "[";
                    if (item.correct && prev_item.correct) {
                        // Correctness came this way
                        out << "color=green";
                    }
                    out << "];" << endl;
                }
            }

        }

        out << "}" << endl;
    }

    out << "}" << endl;
}

Funnel::Item& Funnel::get_item(size_t index) {
    assert(!stages.empty());
    if (index >= stages.back().items.size()) {
        // Allocate up through here
        stages.back().items.resize(index + 1);
    }
    return stages.back().items[index];
}

size_t Funnel::create_item() {
    assert(!stages.empty());
    
    // Work out where to put it
    size_t next_index = stages.back().projected_count;
    // Make sure the item slot exists
    get_item(next_index);
    // Record the item's creation
    stages.back().projected_count++;
    
    // Return the index used
    return next_index;
}
    



}













