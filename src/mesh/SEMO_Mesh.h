#pragma once
class SEMO_Mesh {
private:
	SEMO_Mesh();
	~SEMO_Mesh();
	
private:
static SEMO_Mesh* m_pInst;

public:
	static SEMO_Mesh* GetInst();
	static void DestroyInst();
};