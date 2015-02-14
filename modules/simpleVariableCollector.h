#ifndef simpleVariableCollector_h
#define simpleVariableCollector_h

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <TNtuple.h>
#include <TFile.h>

class simpleVariableCollector
{
  public:
  
    simpleVariableCollector () {} ;
    ~simpleVariableCollector () 
      {
        for (std::map<std::string, std::vector<float> * >::iterator it = m_container.begin () ;
             it != m_container.end () ; ++it)
          delete it->second ;
      }      

    int addVariable (std::string varname)
      {
        if (m_container.find (varname) != m_container.end ())
          return 1 ;
        m_container[varname] = new std::vector<float> () ;
        return 0 ;
      }

    int addVariable2D (std::string varname)
      {
        if (m_container2D.find (varname) != m_container2D.end ())
          return 1 ;
        m_container2D[varname] = new std::vector<std::vector<float> > () ;
        return 0 ;
      }

    int fillContainer (std::string varname, float value)
      {
        if (m_container.find (varname) == m_container.end ())
          return 1 ;
        m_container[varname]->push_back (value) ;
        return 0 ; 
      }

    int fillContanier2D (std::string varname, float first, float second)
      {
        if (m_container2D.find (varname) == m_container2D.end ())
          return 1 ;
        std::vector<float> dummy (2) ;
        m_container2D[varname]->push_back (dummy) ;
        m_container2D[varname]->back ().at (0) = first ;
        m_container2D[varname]->back ().at (1) = second ;
        return 0 ; 
      }

    void save (std::string outputFileName)
      {
        std::cout << "opening " << outputFileName << " for saving..." << std::endl ;
        TFile outfile (outputFileName.c_str (), "recreate") ;
        outfile.cd () ;
        for (std::map<std::string, std::vector<float> * >::iterator it = m_container.begin () ;
             it != m_container.end () ; ++it)
          {
            TString name = it->first.c_str () ;
            TNtuple N (name, name, name) ;
            for (unsigned int i = 0 ; i < it->second->size () ; ++i)
              N.Fill (it->second->at (i)) ;
            N.Write () ;  
          }
        for (std::map<std::string, std::vector<std::vector<float>> * >::iterator it = m_container2D.begin () ;
             it != m_container2D.end () ; ++it)
          {
            TString name = it->first.c_str () ;
            name.ReplaceAll (":", "_") ;
            TNtuple N (name, name, "fir:sec") ;
            for (unsigned int i = 0 ; i < it->second->size () ; ++i)
              N.Fill (it->second->at (i).at (0), it->second->at (i).at (1)) ;
            N.Write () ;  
          }
        outfile.Close () ;
        return ;  
      }
      
      
  private:
  
    std::map<std::string, std::vector<float> * >  m_container ;
    std::map<std::string, std::vector<std::vector<float> > * >  m_container2D ;

} ;


#endif