#ifndef DestructiveCSUI_H
#define DestructiveCSUI_H



#include "DestructiveCSSceneObject.h"
#include "Scene.h"

#include "UIManager.h"

#define SHOW_FORCE_SLIDER(MODE)                                   \
    if (currentAnimation == (MODE)) {                             \
        ImGui::SameLine();                                        \
        ImGui::PushItemWidth(50.0f);                              \
        char sliderLabel[32];                                     \
        snprintf(sliderLabel, sizeof(sliderLabel), "Magnitude##%d", (MODE)); \
        if (ImGui::SliderFloat(sliderLabel, &forceMagnitude, 0.0f, 2.0f)) { \
            systemRef->applyExternalForces((MODE), forceMagnitude); \
        }                                                         \
        ImGui::PopItemWidth();                                    \
    }




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
        for (const auto & sceneObject : SceneManager::GetInstance().getSceneObjects()) {
            DestructiveCSSceneObject * destructiveCSSceneObject = dynamic_cast<DestructiveCSSceneObject*>(sceneObject.get());
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
        const char* sceneNames[] = { "Bunny", "Armadillo", "GummyBear", "Skycraper", "Teapot" };
        ImGui::SetNextWindowPos(ImVec2(5, 10), ImGuiCond_Always);
		ImGui::SetNextWindowSize(ImVec2(340, ImGui::GetIO().DisplaySize.y-20), ImGuiCond_Always);
        if (ImGui::Begin("DestructiveSystem", nullptr, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize)) {
            ImGui::BeginChild("SystemInfo", ImVec2(0,150), true);
			ImGui::Text("Particle Count: %d", systemRef->getNumParticles());
            ImGui::Text("Voxel Count: %d", systemRef->getNumVoxels());
            ImGui::Text("Face Count: %d", systemRef->getNumFaces());
            ImGui::Text("Active Voxel Count: %d", systemRef->getNumVoxels());
            ImGui::Text("Active Face Count: %d", systemRef->getNumFaces());
            ImGui::Text("Num Substeps: %d", systemRef->getSubsteps());
            ImGui::Text("Num Iterations: %d", 1);
            ImGui::EndChild();

            ImGui::BeginChild("SystemGlobal", ImVec2(0,90), true);
            if (ImGui::Button("Reinitialize Scene")) {
                // reinitialize the system
                if(systemRef){
                    systemRef->reinitializeSystem();
                    systemRef->applyExternalForces(currentAnimation, forceMagnitude);
                }
                else{
                    std::cerr << "Error: [DestructiveCSUI_left] systemRef is nullptr" << std::endl;
                }
            }
            if (ImGui::Checkbox("GPU simulation", &bSystemComputeGPU)) {
                systemRef->setComputeComponentActive(bSystemComputeGPU);
            }
            ImGui::SameLine();
            ImGui::PushItemWidth(100.0f); 
            if (ImGui::SliderFloat("DT", &dt, 3.0f, 30.0f)) {
                systemRef->setDT(dt/10000.0f);
            }
            ImGui::PopItemWidth();
            if (ImGui::Combo("Select Scene", &currentScene, sceneNames, IM_ARRAYSIZE(sceneNames))) {
                if (currentScene != previousScene) {
                    previousScene = currentScene;
                    
                    SceneManager::GetInstance().switchScene(static_cast<SceneName>(currentScene));
                    currentAnimation = 0;
                    forceMagnitude = 1.0f;
                    projectileType = 0;
                    bool success = ObtainDestructiveCSSceneObject();
                    if (!success) {
                       std::cerr << "Error: [DestructiveCSUI_left] failed to obtain DestructiveCSSceneObject" << std::endl;
                    }
                    applyAllParams();
                }
            }
            ImGui::EndChild();
            ImGui::BeginChild("SystemVisual", ImVec2(0,110), true);
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
            ImGui::EndChild();

            ImGui::BeginChild("Animation", ImVec2(0,370), true);
            ImGui::Text("Animation Type");
            if (ImGui::RadioButton("None", currentAnimation == 0)) {
                currentAnimation = 0;
                forceMagnitude = 1.0f;
                systemRef->applyExternalForces(0, 0.0f);
            }
            if (ImGui::RadioButton("Explode", currentAnimation == 1)) {
                currentAnimation = 1;
                forceMagnitude = 1.0f;
                systemRef->applyExternalForces(1, 1.0f);
            }
            SHOW_FORCE_SLIDER(1);
            if (ImGui::RadioButton("Spin", currentAnimation == 2)) {
                currentAnimation = 2;
                forceMagnitude = 1.0f;
                systemRef->applyExternalForces(2, 1.0f);
            }
            SHOW_FORCE_SLIDER(2);
            if (ImGui::RadioButton("Split X", currentAnimation == 3)) {
                currentAnimation = 3;
                forceMagnitude = 1.0f;
                systemRef->applyExternalForces(3, 1.0f);
            }
            SHOW_FORCE_SLIDER(3);
            if (ImGui::RadioButton("Split Y", currentAnimation == 4)) {
                currentAnimation = 4;
                forceMagnitude = 1.0f;
                systemRef->applyExternalForces(4, 1.0f);
            }
            SHOW_FORCE_SLIDER(4);
            if (ImGui::RadioButton("Split Z", currentAnimation == 5)) {
                currentAnimation = 5;
                forceMagnitude = 1.0f;
                systemRef->applyExternalForces(5, 1.0f);
            }
            SHOW_FORCE_SLIDER(5);
            if (ImGui::RadioButton("Vortex X", currentAnimation == 6)) {
                currentAnimation = 6;
                forceMagnitude = 1.0f;
                systemRef->applyExternalForces(6, 1.0f);
            }
            SHOW_FORCE_SLIDER(6);
            if (ImGui::RadioButton("Vortex Y", currentAnimation == 7)) {
                currentAnimation = 7;
                forceMagnitude = 1.0f;
                systemRef->applyExternalForces(7, 1.0f);
            }
            SHOW_FORCE_SLIDER(7);
            if (ImGui::RadioButton("Vortex Z", currentAnimation == 8)) {
                currentAnimation = 8;
                forceMagnitude = 1.0f;
                systemRef->applyExternalForces(8, 1.0f);
            }
            SHOW_FORCE_SLIDER(8);
            if (ImGui::RadioButton("Turb", currentAnimation == 9)) {
                currentAnimation = 9;
                forceMagnitude = 1.0f;
                systemRef->applyExternalForces(9, 1.0f);
            }
            SHOW_FORCE_SLIDER(9);
            if (ImGui::RadioButton("Shear", currentAnimation == 10)) {
                currentAnimation = 10;
                forceMagnitude = 1.0f;
                systemRef->applyExternalForces(10, 1.0f);
            }
            SHOW_FORCE_SLIDER(10);
            if (ImGui::RadioButton("Stretch X", currentAnimation == 11)) {
                currentAnimation = 11;
                forceMagnitude = 1.0f;
                systemRef->applyExternalForces(11, 1.0f);
            }
            SHOW_FORCE_SLIDER(11);
            if (ImGui::RadioButton("Stretch Y", currentAnimation == 12)) {
                currentAnimation = 12;
                forceMagnitude = 1.0f;
                systemRef->applyExternalForces(12, 1.0f);
            }
            SHOW_FORCE_SLIDER(12);
            if (ImGui::RadioButton("Twist X", currentAnimation == 13)) {
                currentAnimation = 13;
                forceMagnitude = 1.0f;
                systemRef->applyExternalForces(13, 1.0f);
            }
            SHOW_FORCE_SLIDER(13);
            if (ImGui::RadioButton("Twist Y", currentAnimation == 14)) {
                currentAnimation = 14;
                forceMagnitude = 1.0f;
                systemRef->applyExternalForces(14, 1.0f);
            }
            SHOW_FORCE_SLIDER(14);  
            ImGui::EndChild();
            ImGui::BeginChild("Projectile", ImVec2(0,100), true);
            ImGui::Text("Projectile Type");
            if (ImGui::RadioButton("None", projectileType == 0)) {
                projectileType = 0;
                systemRef->setProjectileType(0);
            }
            if (ImGui::RadioButton("Sphere linear", projectileType == 1)) {
                projectileType = 1;
                systemRef->setProjectileType(1);
            }
            if (ImGui::RadioButton("Box blender", projectileType == 2)) {
                projectileType = 2;
                systemRef->setProjectileType(2);
            }
            ImGui::EndChild();
                
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
    int currentAnimation = 0;
    int projectileType = 0;
    float forceMagnitude = 1.0f;
    float dt = 15.0f;
    DestructiveCSSceneObject * systemRef;

    void applyAllParams() {
        systemRef->setComputeComponentActive(bSystemComputeGPU);
        systemRef->setRenderSkinMesh(bDrawMesh);
        systemRef->setRenderParticles(bDrawParticle);
        systemRef->setRenderVoxels(bDrawVoxel);
        systemRef->setRenderOriginalMesh(bDrawOriginalMesh);
        systemRef->setDT(dt/10000.0f);
    }
};





#endif // DestructiveCSSceneObject_H
