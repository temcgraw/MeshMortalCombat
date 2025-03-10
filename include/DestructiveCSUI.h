#ifndef DestructiveCSUI_H
#define DestructiveCSUI_H



#include "DestructiveCSSceneObject.h"
#include "Scene.h"

#include "UIManager.h"
// the UI for such system, I put it here
class DestructiveCSUI_left : public UIWindow {
public:
    DestructiveCSUI_left() {
        bool success = ObtainDestructiveCSSceneObject();
        if (!success) {
            std::cerr << "Error: [DestructiveCSUI_left] failed to obtain DestructiveCSSceneObject" << std::endl;
        }
        applyAllParams();
    }
    bool ObtainDestructiveCSSceneObject(){
        DestructiveCSSceneObject * destructiveCSSyetem = nullptr;
        // try find the DestructiveCSSceneObject in the scene objects
        for (SceneObject * sceneObject : SceneManager::GetInstance().getSceneObjects()) {
            DestructiveCSSceneObject * destructiveCSSceneObject = dynamic_cast<DestructiveCSSceneObject*>(sceneObject);
            if (destructiveCSSceneObject) {
                destructiveCSSyetem = destructiveCSSceneObject;
                break;
            }
        }
        if (!destructiveCSSyetem) {
            return false;
        }
        systemRef = destructiveCSSyetem;
        return true;
    }

    void drawWindow() override {
        if (!display) return;
        const char* sceneNames[] = { "SCENE1", "SCENE2", "SCENE3", "SCENE4", "SCENE5" };
        ImGui::SetNextWindowPos(ImVec2(5, 10), ImGuiCond_Always);
		ImGui::SetNextWindowSize(ImVec2(340, ImGui::GetIO().DisplaySize.y-20), ImGuiCond_Always);
        if (ImGui::Begin("DestructiveSystem", nullptr, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize)) {
			ImGui::Text("Particle Count: %d", systemRef->getNumParticles());
            ImGui::Text("Voxel Count: %d", systemRef->getNumVoxels());
            ImGui::Text("Face Count: %d", systemRef->getNumFaces());
            ImGui::Text("Active Voxel Count: %d", systemRef->getNumVoxels());
            ImGui::Text("Active Face Count: %d", systemRef->getNumFaces());
            ImGui::Text("Num Substeps: %d", systemRef->getSubsteps());
            ImGui::Text("Num Iterations: %d", 1);

            if (ImGui::Button("Reinitialize Scene")) {
                // reinitialize the system
                if(systemRef){
                    systemRef->reinitializeSystem();
                }
                else{
                    std::cerr << "Error: [DestructiveCSUI_left] systemRef is nullptr" << std::endl;
                }
            }
            if (ImGui::Checkbox("GPU simulation", &bSystemComputeGPU)) {
                systemRef->setComputeComponentActive(bSystemComputeGPU);
            }
            if (ImGui::Combo("Scene Selection", &currentScene, sceneNames, IM_ARRAYSIZE(sceneNames))) {
                if (currentScene != previousScene) {
                    previousScene = currentScene;
                    
                    SceneManager::GetInstance().switchScene(static_cast<SceneName>(currentScene));
                        
                    bool success = ObtainDestructiveCSSceneObject();
                    if (!success) {
                       std::cerr << "Error: [DestructiveCSUI_left] failed to obtain DestructiveCSSceneObject" << std::endl;
                    }
                    applyAllParams();
                }
            }
            if (ImGui::Checkbox("Draw Mesh", &bDrawMesh)) {
                systemRef->setRenderSkinMesh(bDrawMesh);
            }
            if (ImGui::Checkbox("Draw Particle", &bDrawParticle)) {
                systemRef->setRenderParticles(bDrawParticle);
            }
            if (ImGui::Checkbox("Draw Voxel", &bDrawVoxel)) {
                systemRef->setRenderVoxels(bDrawVoxel);
            }
            if (ImGui::Checkbox("Draw Original Mesh", &bDrawOriginalMesh)) {
                systemRef->setRenderOriginalMesh(bDrawOriginalMesh);
            }
            ImGui::End();
		}
    }

private:
    bool bSystemComputeGPU = false;
    bool bDrawMesh = true;
    bool bDrawParticle = false;
    bool bDrawVoxel = true;
    bool bDrawOriginalMesh = false;
    int currentScene = 0;
    int previousScene = 0;
    DestructiveCSSceneObject * systemRef;

    void applyAllParams() {
        systemRef->setComputeComponentActive(bSystemComputeGPU);
        systemRef->setRenderSkinMesh(bDrawMesh);
        systemRef->setRenderParticles(bDrawParticle);
        systemRef->setRenderVoxels(bDrawVoxel);
        systemRef->setRenderOriginalMesh(bDrawOriginalMesh);
    }
};





#endif // DestructiveCSSceneObject_H
