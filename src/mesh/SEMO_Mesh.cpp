#include "SEMO_Mesh.h"

SEMO_Mesh* SEMO_Mesh::m_pInst = nullptr;

SEMO_Mesh::SEMO_Mesh(){
}

SEMO_Mesh::~SEMO_Mesh(){
}

SEMO_Mesh* SEMO_Mesh::GetInst(){
	if(!m_pInst){
		m_pInst = new SEMO_Mesh;
	}
	return m_pInst;
}

void SEMO_Mesh::DestroyInst(){
	if(!m_pInst){
		return;
	}
	delete m_pInst;
	m_pInst = nullptr;
}