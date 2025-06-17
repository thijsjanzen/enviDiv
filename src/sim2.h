#include <vector>
#include <string>
#include <array>
#include "random_thijs.h"
#include "util.h"

#include "Rcpp.h"

enum info {btime, parent, id, death};
enum pars {extinction, sym_high, sym_low, allo, wobble, water};
enum ww   {low, high};

class entry {
public:
  entry(double birth, int p, int ID) :
  b_time_(birth),
  parent_id_(p),
  id_(ID){
    d_time_ = -1.0;
    in_num_pockets = 1;
    dead_ = false;
  }

  void die(const double& t) {
    d_time_ = t;
    dead_ = true;
  }

  bool dead() const {
    return dead_;
  }

  bool alive() const {
    return !dead_;
  }

  int parent() const {
    return parent_id_;
  }

  int ID() const {
    return id_;
  }

  bool die_pocket() {
    in_num_pockets--;
    if (in_num_pockets == 0) return true;
    return false;
  }

  double get_btime() const {
    return b_time_;
  }
  double get_dtime() const {
    return d_time_;
  }



  int in_num_pockets;

private:
  double d_time_;
  bool dead_;


  const double b_time_;
  const int parent_id_;
  const int id_;
};

namespace new_sim {

struct simulation {

  double t;
  
  const std::array<double, 4> params_;
  const double max_time;
  const int max_species;
  const int focal_model;

  std::vector<entry> ltable;
  std::array<double, 4> rates;
  ww waterlevel;
  rnd_t rnd;
  std::string run_info;
  int waterlevelchanges;
  double last_w_change;

  std::array<int, 2> crowns;


  simulation(const std::array<double, 4>& p,
             int chosen_model,
             double crown_age,
             int max_num_spec,
             int seed) :
    params_(p),
    max_time(crown_age),
    max_species(max_num_spec),
    focal_model(chosen_model) {
      seed > 0 ? rnd.set_seed(seed) : rnd.set_seed(time(NULL));
  }

  void run() {
    run_info = "not_run_yet";
    t = 0.0;
    // water levels start at t = 0.0
    const std::vector<double> waterlevels = get_waterlevel_changes(focal_model,
                                                                   max_time,
                                                                   rnd,
                                                                   params_[ pars::water ]);

    waterlevelchanges = 1;
    last_w_change = -1;
    auto next_w_change = waterlevels[waterlevelchanges];
    ltable.clear();
    ltable.emplace_back(entry(0.0, 0, -1));
    ltable.emplace_back(entry(0.0, -1, 2));

    waterlevel = high;

    crowns = {1, 1};

    while( true ) {
      update_rates();

      double dt = draw_dt();

     // std::cerr << t << " " << dt << "\n";

      t += dt;

      if (t >= next_w_change) {
        change_water_level(next_w_change);
        t = next_w_change;
        last_w_change = next_w_change;
        waterlevelchanges++;
        next_w_change = waterlevels[waterlevelchanges];
        continue;
      }

      if (t >= max_time)  {
        run_info = "done";
        break;
      }

      pars event = draw_event();
      std::string txt = "place_holder";
      switch(event) {
        case allo: txt =  "allo"; break;
        case sym_high: txt = "sym_high"; break;
        case sym_low: txt = "sym_low"; break;
        case extinction: txt = "extinction"; break;
        case water: txt = "water"; break;
        case wobble: txt = "wobble"; break;
      }
      //std::cerr << txt << "\n";


      apply_event(event);

      if (crowns[0] < 1 || crowns[1] < 1) {
        run_info = "extinct";
        break;
      }

      if (crowns[0] + crowns[1] > max_species) {
        run_info = "overshoot";
        break;
      }
    }
  }

  void apply_event(pars event) {
     switch(event) {
       case extinction: event_extinction(); break;
       case sym_high  : event_sym_high();   break;
       case sym_low   : event_sym_low(last_w_change);    break;
       case allo      : event_allo(last_w_change)   ;    break;
     }
  }

  void death(size_t index) {
    ltable[index].die(t);
    if (ltable[index].ID() < 0) {
      crowns[0]--;
    } else {
      crowns[1]--;
    }
  }

  void birth(double t, double parent_index) {
    auto parent_id = ltable[parent_index].ID();
    int new_id = static_cast<int>(ltable.size()) + 1; // start counting at 1, add one
    if (parent_id < 0) {
      new_id *= -1;
      crowns[0]++;
    } else {
      crowns[1]++;
    }

    ltable.emplace_back(entry(t, parent_id, new_id));
  }


  void event_extinction() {
    size_t index = draw_extinct_species();

    if (waterlevel == high) {
      death(index);
    } else {
      bool local_extinction = ltable[index].die_pocket();
      if (local_extinction) death(index);
    }
  }

  void event_sym_high() {
    size_t index = rnd.random_number(ltable.size());
    while(ltable[index].dead()) index = rnd.random_number(ltable.size());

    birth(t, index);
  }

  void event_sym_low(const double& last_waterlevel_change) {
    size_t index = rnd.random_number(ltable.size());
    while(ltable[index].dead()) index = rnd.random_number(ltable.size());

    if (ltable[index].in_num_pockets == 1) {
      // 'normal' sympatric speciation
      birth(t, index);
    } else {
      // sympatric speciation in one pocket, but not in the other (!)
      // remove parent species from one pocket
      ltable[index].in_num_pockets = 1;

      birth(last_waterlevel_change, index);

      birth(t, ltable.size() - 1);
    }
  }

  void event_allo(const double& last_waterlevel_change) {
    size_t index = rnd.random_number(ltable.size());
    while(ltable[index].dead() && ltable[index].in_num_pockets != 2) {
      index = rnd.random_number(ltable.size());
    }

    ltable[index].in_num_pockets = 1;

    birth(last_waterlevel_change, index);
  }

  void update_rates() {
    rates = {0.0, 0.0, 0.0, 0.0};
    for (const auto& i : ltable) {
      if (i.alive()) {
        rates[extinction] += params_[extinction] * i.in_num_pockets;
        if (waterlevel == high) {
          rates[sym_high] += params_[sym_high];
        } else {
          rates[sym_low] += params_[sym_low] * i.in_num_pockets;

          if (i.in_num_pockets == 2) {
            rates[allo]    += params_[allo];
          }
        }
      }
    }
  }

  pars draw_event() {
    double total_rate = std::accumulate(rates.begin(), rates.end(), 0.0);
    double r = rnd.uniform(0.0, total_rate);
    for (size_t i = 0; i < rates.size(); ++i) {
      r -= rates[i];
      if (r <= 0.0) {
        return static_cast<pars>(i);
      }
    }
    return pars::allo;
  }

  double draw_dt() {
    double total_rate = std::accumulate(rates.begin(), rates.end(), 0.0);
    return rnd.Expon(total_rate);
  }

  void change_water_level(double w_t) {
    if (waterlevel == high) {
      waterlevel = low;
    } else {
      waterlevel = high;
    }

    for (auto& i : ltable) {
      if (i.alive()) {
        if (waterlevel == low) {
          i.in_num_pockets = 2;
        } else {
          i.in_num_pockets = 1;
        }
      }
    }
  }

  size_t draw_extinct_species() {
    size_t index = rnd.random_number(ltable.size());
    while( true ) {
      if (ltable[index].dead()) continue;

      double prob = ltable[index].in_num_pockets * 0.5;
      if (rnd.uniform() < prob) break;

      index = rnd.random_number(ltable.size());
    }
    return index;
  }

  Rcpp::NumericMatrix get_ltable() {
    Rcpp::NumericMatrix out(ltable.size(), 4);
    for (size_t i = 0; i < ltable.size(); ++i) {
      out(i, 0) = ltable[i].get_btime();
      out(i, 1) = ltable[i].parent();
      out(i, 2) = ltable[i].ID();
      out(i, 3) = ltable[i].get_dtime();
    }
    return out;
  }
};

} // ned namespace new_sim
