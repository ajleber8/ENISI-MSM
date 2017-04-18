#include "NeutrophilGroup.h"
#include "agent/ENISIAgent.h"
#include "compartment/Compartment.h"
#include "grid/Properties.h"
#include "DataWriter/LocalFile.h"
#include <iostream>
#include <iomanip>

using namespace ENISI;

NeutrophilGroup::NeutrophilGroup(Compartment * pCompartment, const double & concentrations) :
										  GroupInterface(pCompartment)
{
	size_t LocalCount = mpCompartment->localCount(concentrations);

	for (size_t i = 0; i < LocalCount; i++){
		mpCompartment->addAgentToRandomLocation(new Agent(Agent::Neutrophil, NeutrophilState::BASE));
	}

	const Properties * pModel = Properties::instance(Properties::model);

	pModel->getValue("p_nrecmac", p_nrecmac);
	pModel->getValue("p_nrecei", p_nrecei);
	pModel->getValue("p_nreced", p_nreced);
	pModel->getValue("p_nbactact", p_nbactact);
	pModel->getValue("p_nepiact", p_nepiact);
	pModel->getValue("p_nbasedeath", p_nbasedeath);
	pModel->getValue("p_nskipecoli", p_nskipecoli);
	pModel->getValue("p_nkillbac", p_nkillbac);
	pModel->getValue("p_nactdeath", p_nactdeath);
	pModel->getValue("p_th1max", p_th1max);
	pModel->getValue("p_neutbaserec", p_neutbaserec);
	pModel->getValue("p_Nbasal", p_Nbasal);										  
											 
}

void NeutrophilGroup::act(const repast::Point<int> & pt)
{
	std::vector< double > Location(2, 0);

	std::vector< Agent * > Neutrophils;
	mpCompartment->getAgents(pt, Agent::Neutrophil, Neutrophils);

	std::vector< Agent * > BacteriaDAs;
	std::vector< Agent * > Tcells;
	std::vector< Agent * > EpithelialCells;
	std::vector< Agent * > Macrophages;

	if (mpCompartment->getType() == Compartment::lamina_propria)
	{
		mpCompartment->getAgents(pt, Agent::BacteriaDA, BacteriaDAs);
		mpCompartment->getAgents(pt, Agent::Macrophage, Macrophages);
		mpCompartment->getAgents(pt, Agent::Tcell, Tcells);
	}
	if (mpCompartment->gridBorders()->distanceFromBorder(pt.coords(), Borders::Y, Borders::LOW) < 0.5)
	{
		mpCompartment->getAgents(pt, 0, -1, Agent::EpithelialCell, EpithelialCells);
	}

	Concentration NeutrophilConcentration;
	concentrations(Agent::Neutrophil, Neutrophils, NeutrophilConcentration);
	Concentration BacteriaDAConcentration;
	concentrations(Agent::BacteriaDA, BacteriaDAs, BacteriaDAConcentration);
	Concentration MacrophageConcentration;
	concentrations(Agent::Macrophage, Macrophages, MacrophageConcentration);
	Concentration TcellConcentration;
	concentrations(Agent::Tcell, Tcells, TcellConcentration);
	Concentration EpithelialCellConcentration;
	concentrations(Agent::EpithelialCell, EpithelialCells, EpithelialCellConcentration);

	double neutConcentration = NeutrophilConcentration[NeutrophilState::BASE] + NeutrophilConcentration[NeutrophilState::ACTIVATED];
	double bacteriaDAConcentration = BacteriaDAConcentration[BacteriaDAState::ECOLI] + BacteriaDAConcentration[BacteriaDAState::MYCO] + BacteriaDAConcentration[BacteriaDAState::KLEB] + BacteriaDAConcentration[BacteriaDAState::ENTERO] + BacteriaDAConcentration[BacteriaDAState::CORIO];
	double ecoliConcentration = BacteriaDAConcentration[BacteriaDAState::ECOLI];
	double infmacConcentration = MacrophageConcentration[MacrophageState::INFLAMMATORY];
	double th1Concentration = TcellConcentration[TcellState::TH1];
	double th17Concentration = TcellConcentration[TcellState::TH17];
	double damagedEpithelialCellConcentration = EpithelialCellConcentration[EpithelialCellState::DAMAGED];
	double inflamedEpithelialCellConcentration = EpithelialCellConcentration[EpithelialCellState::INFLAMED];

	std::vector< Agent * >::iterator it = Neutrophils.begin();
	std::vector< Agent * >::iterator end = Neutrophils.end();

	for (; it != end; ++it)
	{
		Agent * pAgent = *it;
		NeutrophilState::State state = (NeutrophilState::State) pAgent->getState();
		NeutrophilState::State newState = state;
		
		if (state == NeutrophilState::DEAD) continue;

		if (state == NeutrophilState::BASE)
		{
			// neutrophils recruited by macrophages
			if (th17Concentration > ENISI::Threshold
					&& p_nrecmac > repast::Random::instance()->createUniDoubleGenerator(0.0, 1.0).next())
			{
				mpCompartment->getLocation(pAgent->getId(), Location);
				mpCompartment->addAgent(new Agent(Agent::Neutrophil, pAgent->getState()), Location);
			}
			// neutrophils recruited by damaged epithelial cells
			if (inflamedEpithelialCellConcentration > ENISI::Threshold
					&& p_nrecei > repast::Random::instance()->createUniDoubleGenerator(0.0, 1.0).next())
			{
				mpCompartment->getLocation(pAgent->getId(), Location);
				mpCompartment->addAgent(new Agent(Agent::Neutrophil, pAgent->getState()), Location);
			}
			// neutrophils recruited by damaged epithelial cells
			if (damagedEpithelialCellConcentration > ENISI::Threshold
					&& p_nreced > repast::Random::instance()->createUniDoubleGenerator(0.0, 1.0).next())
			{
				mpCompartment->getLocation(pAgent->getId(), Location);
				mpCompartment->addAgent(new Agent(Agent::Neutrophil, pAgent->getState()), Location);
			}
			/*if (p_neutbaserec > repast::Random::instance()->createUniDoubleGenerator(0.0, 1.0).next())
			{
				mpCompartment->getLocation(pAgent->getId(), Location);
				mpCompartment->addAgent(new Agent(Agent::Neutrophil, pAgent->getState()), Location);
			}*/
			if (p_Nbasal > neutConcentration)
			{
				mpCompartment->getLocation(pAgent->getId(), Location);
				mpCompartment->addAgent(new Agent(Agent::Neutrophil, pAgent->getState()), Location);
			}
			// neutrophils activated by bacteria
			if (bacteriaDAConcentration > ENISI::Threshold
					&& p_nbactact > repast::Random::instance()->createUniDoubleGenerator(0.0,1.0).next())
			{
				newState = NeutrophilState::ACTIVATED;
				pAgent->setState(newState);
			}
			// neutrophils activated by epithelial cells
			if (inflamedEpithelialCellConcentration > ENISI::Threshold
					&& p_nepiact > repast::Random::instance()->createUniDoubleGenerator(0.0,1.0).next())
			{
				newState = NeutrophilState::ACTIVATED;
				pAgent->setState(newState);
			}
			// neutrophils die
			if (p_nbasedeath > repast::Random::instance()->createUniDoubleGenerator(0.0,1.0).next())
			{
				mpCompartment->removeAgent(pAgent);
			}
		}

		if (state == NeutrophilState::ACTIVATED)
		{
			// activated neutrophil efficacy is diminished by ecoli
			if (ecoliConcentration > ENISI::Threshold
					&& p_nskipecoli > repast::Random::instance()->createUniDoubleGenerator(0.0,1.0).next())
			{
				continue;
			}
			// activated neutrophils kill bacteria and secrete cytokines
			else if (bacteriaDAConcentration > ENISI::Threshold
					&& p_nkillbac > repast::Random::instance()->createUniDoubleGenerator(0.0,1.0).next())
			{
				mpCompartment->cytokineValue("eIL17", pt) += 1;
				if (BacteriaDAs.size() > 0)
				{
					mpCompartment->removeAgent(BacteriaDAs[BacteriaDAs.size() - 1]);
					BacteriaDAs.pop_back();
				}
			}
			// activated neutrophils die
			if (p_nactdeath > repast::Random::instance()->createUniDoubleGenerator(0.0,1.0).next() && th1Concentration * repast::Random::instance()->createUniDoubleGenerator(0.0,1.0).next() < p_th1max)
			{
				mpCompartment->removeAgent(pAgent);
			}
		}
	}
}

// virtual
void NeutrophilGroup::move()
{
	// TODO CRITICAL Determine the maximum speed
	double MaxSpeed = 4.0;
	bool directed [9];
	int max = 0;
	int dir = 0;
	// Find all local agents and move them
	Compartment::LocalIterator itLocal = mpCompartment->localBegin();
	Compartment::LocalIterator endLocal = mpCompartment->localEnd();

	for (; itLocal != endLocal; ++itLocal)
	{
        std::vector< double > Location(2);
        mpCompartment->getLocation((*itLocal)->getId(),Location);
        const repast::Point<int> & pt = mpCompartment->spaceToGrid(Location);

		int count = 0;
		for (int x = -1; x < 2; ++x) {
			for (int y = -1; y < 2; ++y) {
				if(mpCompartment->cytokineValue("eIL17", pt, x, y) >= max){
					directed[count] = true;
					max = mpCompartment->cytokineValue("eIL17", pt, x, y);
				}
				count++;
			}
			count ++;
		}
		for (int i = 0; i < 9; ++i) {
			if(directed[i]){
				dir = i;
			}
		}
		if(dir == 4){
			mpCompartment->moveRandom((*itLocal)->getId(), MaxSpeed);
		} else{
			mpCompartment->moveRandom((*itLocal)->getId(), MaxSpeed * 0.5);
			mpCompartment->moveDirected((*itLocal)->getId(), MaxSpeed * 0.5,dir);
		}
	}
}

// virtual
void NeutrophilGroup::write(const repast::Point<int> &)
{

}
