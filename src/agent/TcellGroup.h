#ifndef ENISI_MSM_AGENT_TCELL_COMPONENT_H
#define ENISI_MSM_AGENT_TCELL_COMPONENT_H

#include "TcellODE.h"
#include "CellGroup.h"
#include "Cytokines.h"

class TcellGroup: public CellGroup 
{
public:
  enum State { NAIVE, TH1, TH17, TREG, DEAD, LAST_STATE_DO_NOT_MOVE};

  struct StateCount { 
    unsigned int state[LAST_STATE_DO_NOT_MOVE]; 
    StateCount() : state() { }
  };

  typedef std::map<repast::Point<int>, StateCount> CoordMap;
  typedef std::map<int, std::vector<std::pair<int, int> > > Transfers;

  TcellGroup(const boost::uintmax_t, CellLayer * p_layer);

  virtual void act();
  int count();
  StateCount countByState();

  void addStateAt(State, const repast::Point<int> &);

  void transferStateTo(State, const repast::Point<int> &, unsigned int = 1);

  const StateCount &
    getCellsAt(const repast::Point<int> & loc) { return coordMap[loc]; }

  static std::vector<TcellGroup *> & instances()
  {
    static std::vector<TcellGroup *> instances;
    return instances;
  }

  virtual Color getColor() { return black; }
  virtual std::string classname() { return "TcellGroup"; }
protected:
private:
  void act(State, const repast::Point<int> &);
  void init(const boost::uintmax_t);

  CoordMap coordMap;
};

#endif
