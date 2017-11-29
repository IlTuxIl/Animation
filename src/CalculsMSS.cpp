/*
 * CalculsMSS.cpp :
 * Copyright (C) 2016 Florence Zara, LIRIS
 *               florence.zara@liris.univ-lyon1.fr
 *               http://liris.cnrs.fr/florence.zara/
 *
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/** \file CalculsMSS.cpp
Programme calculant pour chaque particule i d un MSS son etat au pas de temps suivant 
 (methode d 'Euler semi-implicite) : principales fonctions de calculs.
\brief Fonctions de calculs de la methode semi-implicite sur un systeme masses-ressorts.
*/ 

#include <stdio.h>
#include <math.h>
#include <vector>
#include <iostream>

#include "vec.h"
#include "ObjetSimule.h"
#include "ObjetSimuleMSS.h"
#include "Viewer.h"

using namespace std;

float calculeNorme(const Vector& a){
    float res = sqrt(pow(a.x, 2) + pow(a.y, 2) + pow(a.z, 2));
    return res;
}

/**
* Calcul des forces appliquees sur les particules du systeme masses-ressorts.
 */
void ObjetSimuleMSS::CalculForceSpring()
{

    Particule* cur;
    Particule* tmp;

    for(int i = 0; i < (int)_SytemeMasseRessort->GetNbParticule(); i++){
        Vector somme;
        cur = _SytemeMasseRessort->GetParticule(i);
        for(int j = 0;  j < (int)_SytemeMasseRessort->_ParticuleList[i]->_RessortList.size(); j++){
            tmp = _SytemeMasseRessort->_ParticuleList[i]->_RessortList[j]->GetParticuleA();
            if(tmp == cur){
                tmp = _SytemeMasseRessort->_ParticuleList[i]->_RessortList[j]->GetParticuleB();
            }

            float raideur = _SytemeMasseRessort->_ParticuleList[i]->_RessortList[j]->GetRaideur();
            float lenRepos = _SytemeMasseRessort->_ParticuleList[i]->_RessortList[j]->GetLrepos();
            float amorti = _SytemeMasseRessort->_ParticuleList[i]->_RessortList[j]->GetAmortissement();

            Vector eij;
            Vector vij;
            eij = -raideur * (calculeNorme(P[cur->_Id] - P[tmp->_Id]) - lenRepos) * (P[cur->_Id] - P[tmp->_Id]) / calculeNorme(P[cur->_Id] - P[tmp->_Id]);
            if(calculeNorme(P[cur->_Id] - P[tmp->_Id]) > 5*lenRepos){
                _SytemeMasseRessort->_ParticuleList[i]->_RessortList.erase(_SytemeMasseRessort->_ParticuleList[i]->_RessortList.begin() + j);
            }
            vij = amorti * cross((V[cur->_Id] - V[tmp->_Id]), (P[cur->_Id] - P[tmp->_Id]) / calculeNorme(P[cur->_Id] - P[tmp->_Id]));
            somme = somme + eij + vij;
        }
        Force[cur->_Id] = somme;
    }
	/// f = somme_i (ki * (l(i,j)-l_0(i,j)) * uij ) + (nuij * (vi - vj) * uij) + (m*g) + force_ext
	
	/// Rq : Les forces dues a la gravite et au vent sont ajoutees lors du calcul de l acceleration
    
		
}//void


/**
 * Gestion des collisions avec le sol - plan (x,y,z).
 */
void ObjetSimuleMSS::CollisionPlan(float x, float y, float z)
{
    /// Arret de la vitesse quand touche le plan
    for(int i = 0; i < P.size(); i++){
        if(P[i].y <= y){
            if(A[i].y < 0){
                A[i].y = 0;
                V[i].y = 0;
            }
        }
    }
    
}// void

