#ifndef UI_MANAGER_H
#define UI_MANAGER_H

#include <Renderer.h>
#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"
#include <FrameRateMonitor.h>
#include <vector>
#include <functional>
#include <string>



class UIWindow {
public:
    std::string name;
    virtual void drawWindow() = 0;
    bool display = true;
};
// UIManager, a class for managing the user interface
// This class is responsible for initializing the ImGui library, rendering the registered windows, and registering new windows
class UIManager {
public:
    void Initialize(GLFWwindow* window){
        this->glfw_window = window;
         //imgui config----------------------
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        // set color theme
        ImGui::StyleColorsLight();
        // embed imgui into glfw and opengl3
        ImGui_ImplGlfw_InitForOpenGL(window, true);
        // OpenGL version 4.3
        ImGui_ImplOpenGL3_Init("#version 430");
        ImGuiIO& io = ImGui::GetIO();
        (void)io;
        // font setting
        io.Fonts->AddFontFromFileTTF("resource/fonts/Cousine-Regular.ttf", 13.0f, NULL, io.Fonts->GetGlyphRangesDefault());
        io.Fonts->AddFontFromFileTTF("resource/fonts/DroidSans.ttf", 13.0f, NULL, io.Fonts->GetGlyphRangesDefault());
        io.Fonts->AddFontFromFileTTF("resource/fonts/Karla-Regular.ttf", 13.0f, NULL, io.Fonts->GetGlyphRangesDefault());
        io.Fonts->AddFontFromFileTTF("resource/fonts/ProggyClean.ttf", 13.0f, NULL, io.Fonts->GetGlyphRangesDefault());
        io.Fonts->AddFontFromFileTTF("resource/fonts/Roboto-Medium.ttf", 13.0f, NULL, io.Fonts->GetGlyphRangesDefault());
    }
    void Render(){
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        for (UIWindow * uiWindow  : windows) {
            uiWindow->drawWindow();
        }

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    }

    void enableAllWindows(){
        for (UIWindow * uiWindow  : windows) {
            uiWindow->display = true;
        }
    }

    void disableAllWindows(){
        for (UIWindow * uiWindow  : windows) {
            uiWindow->display = false;
        }
    }

    // register a window with a name and a draw function
    void RegisterWindow(UIWindow * window){
        windows.push_back(window);
    }
    

private:

    std::vector<UIWindow*> windows; // store all registered windows
    GLFWwindow* glfw_window = nullptr;
};


// --------------------------------------------------------------------------------
// some windows that can be registered for current project
class LogWindow : public UIWindow {
public:
    LogWindow(FrameRateMonitor * _frameRateMonitor, Camera * _camera, bool * _enableSceneTick) {
        frameRateMonitor = _frameRateMonitor;
        camera = _camera;
        enableSceneTick = _enableSceneTick;
        GPU_name = glGetString(GL_RENDERER);
        
    }
    void drawWindow() override {
        if (!display) return;
        ImGui::SetNextWindowPos(ImVec2(ImGui::GetIO().DisplaySize.x - 345, 10), ImGuiCond_Always);
		ImGui::SetNextWindowSize(ImVec2(340, 250), ImGuiCond_Always);
        if (ImGui::Begin("LOG", nullptr, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize)) {

			ImGui::Text("Average FPS: %.1f ", frameRateMonitor->getFPS()); ImGui::SameLine(); ImGui::Text("| %.1f ", frameRateMonitor->getAverageFPS());
            ImGui::Text("Frame Time: %.6f ms", 999);
			ImGui::Text("Camera POS: %.3f %.3f %.3f", camera->Position[0], camera->Position[1], camera->Position[2]);
			ImGui::Text("Camera DIR: %.3f %.3f %.3f", camera->Front[0], camera->Front[1], camera->Front[2]);
			ImGui::Text("Camera FOV: %.3f", camera->Zoom);
            ImGui::Text("Particle Count: %d", 999);
            ImGui::Text("Voxel Count: %d", 999);
            ImGui::Text("Voxel LOD0 Count: %d", 9999);
            ImGui::Text("Num Substeps: %d", 999);
            ImGui::Text("Num Iterations: %d", 999);
            ImGui::Text("Device: %s", GPU_name);
            if (*enableSceneTick) {
                if (ImGui::Button("Disable Scene Tick")) {
                    *enableSceneTick = false;
                }
            } 
            else {
                if (ImGui::Button("Enable Scene Tick")) {
                    *enableSceneTick = true;
                }
            }
            ImGui::End();
		}
    }
private:
    FrameRateMonitor * frameRateMonitor;
    Camera * camera;
    bool * enableSceneTick = nullptr;
    const char * enableSceneTickText = "Enable Scene Tick";
    const char * disableSceneTickText = "Disable Scene Tick";
    const GLubyte* GPU_name = nullptr;
};









#endif