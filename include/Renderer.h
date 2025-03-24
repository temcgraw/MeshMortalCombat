#ifndef RENDERER_H
#define RENDERER_H

#include <glm/glm.hpp>
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <Shader.h>
#include <Camera.h>
#include <Scene.h>
#include <FrameRateMonitor.h>


#include <iostream>
#include <chrono>
#include <list>





struct cameraUBOData {
    glm::mat4 view;
    glm::mat4 projection;
    glm::vec3 cameraPos;
};

struct alignas(16) lightUBOData {
    alignas(4) bool bUsePointLight;
    alignas(4) int PointLightCount;
    alignas(16) glm::vec3 PointLightPos1;
    alignas(16) glm::vec3 PointLightColor1;
    alignas(16) glm::vec3 PointLightPos2;
    alignas(16) glm::vec3 PointLightColor2;
    alignas(16) glm::vec3 PointLightPos3;
    alignas(16) glm::vec3 PointLightColor3;
    alignas(4) bool bUseDirectionalLight;
    alignas(16) glm::vec3 DirectionalLightDir;
    alignas(16) glm::vec3 DirectionalLightColor;
    alignas(16) glm::vec3 AmbientLightColor;
};


class Renderer {
    public:

    GLFWwindow * window = nullptr;
    Shader * shader = nullptr;
    Camera * camera = nullptr;
    int screen_width = 800;
    int screen_height = 600;
    FrameRateMonitor * frameRateMonitor = nullptr;

    GLuint camera_ubo, light_ubo;
    cameraUBOData CameraUploadData;
    lightUBOData LightUploadData;


	Renderer(int _screen_width, int _screen_height, Camera * camera = nullptr) : camera(camera), screen_width(_screen_width), screen_height(_screen_height) {
		initialize();
	}

	void initialize() {
        //----------------------------------------------------------------------------------
        // glfw: initialize and configure
        // ------------------------------
        glfwInit();
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

        glfwWindowHint(GLFW_DOUBLEBUFFER, GLFW_TRUE); // Enable double buffering
        glfwWindowHint(GLFW_SAMPLES, 4); // Enable multisampling



        #ifdef __APPLE__
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
        #endif

        // glfw window creation
        // --------------------
        std::cout << "Creating GLFW window" << std::endl;
        window = glfwCreateWindow(screen_width, screen_height, "DestructiveAdaptiveGrid", NULL, NULL);


        if (window == NULL)
        {
            std::cout << "Failed to create GLFW window" << std::endl;
            glfwTerminate();
            return;
        }
        // set window position and size
        glfwMakeContextCurrent(window);

        glfwSwapInterval(_vSync ? 1 : 0); // Enable vsync

        // tell GLFW to capture our mouse
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);

        // glad: load all OpenGL function pointers
        // ---------------------------------------
        if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
        {
            std::cout << "Failed to initialize GLAD" << std::endl;
            return;
        }


        // configure global opengl state (they are not global, but we treat them as global since we probably won't change them)
        // -----------------------------
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
        glPointSize(8.0);
        glLineWidth(4.0);
        glEnable(GL_MULTISAMPLE);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);

        InitializeUbo();
        // set up the screenCanvas, it will be used to display the ray tracing result (either CPU version or GPU version)


	}

    void InitializeUbo() {
        
        glGenBuffers(1, &camera_ubo);
        glBindBuffer(GL_UNIFORM_BUFFER, camera_ubo);

        CameraUploadData.cameraPos = camera->Position;
        CameraUploadData.projection = glm::perspective(glm::radians(camera->Zoom), static_cast<float>(screen_width) / static_cast<float>(screen_height), 0.1f, 500.0f);
        CameraUploadData.view = camera->GetViewMatrix();
        glBufferData(GL_UNIFORM_BUFFER, sizeof(CameraUploadData), &CameraUploadData, GL_STATIC_DRAW);
        glBindBuffer(GL_UNIFORM_BUFFER, 0);

        GLuint cameraBindingPoint = 0;
        glBindBufferBase(GL_UNIFORM_BUFFER, cameraBindingPoint, camera_ubo);

        LightUploadData.bUsePointLight = false;
        LightUploadData.PointLightCount = 1;
        LightUploadData.PointLightPos1 = glm::vec3(0.0f, 1.0f, 0.0f);
        LightUploadData.PointLightColor1 = glm::vec3(0.3f, 0.5f, 0.3f);
        LightUploadData.bUseDirectionalLight = true;
        LightUploadData.DirectionalLightDir = glm::vec3(1.0f, -1.0f, -0.5f);
        LightUploadData.DirectionalLightColor = glm::vec3(0.9f, 0.9f, 1.0f);
        LightUploadData.AmbientLightColor = glm::vec3(0.4,0.4,0.5);

        glGenBuffers(1, &light_ubo);
        glBindBuffer(GL_UNIFORM_BUFFER, light_ubo);
        glBufferData(GL_UNIFORM_BUFFER, sizeof(LightUploadData), &LightUploadData, GL_STATIC_DRAW);
        glBindBuffer(GL_UNIFORM_BUFFER, 0);

        GLuint lightBindingPoint = 10;
        glBindBufferBase(GL_UNIFORM_BUFFER, lightBindingPoint, light_ubo);


    }

    void uploadUbo() {
        if (camera_ubo == 0) {
            glGenBuffers(1, &camera_ubo);
            glBindBuffer(GL_UNIFORM_BUFFER, camera_ubo);
        }
        else {
            glBindBuffer(GL_UNIFORM_BUFFER, camera_ubo);
        }
        glBufferData(GL_UNIFORM_BUFFER, sizeof(CameraUploadData), &CameraUploadData, GL_STATIC_DRAW);
        glBindBuffer(GL_UNIFORM_BUFFER, 0);
    }

    void updateUboData() {
		CameraUploadData.cameraPos = camera->Position;
		CameraUploadData.projection = glm::perspective(glm::radians(camera->Zoom), static_cast<float>(screen_width) / static_cast<float>(screen_height), 0.1f, 500.0f);
		CameraUploadData.view = camera->GetViewMatrix();
	}

    void poll_events() {
        glfwPollEvents();
    }

    void swap_buffers() {
		glfwSwapBuffers(window);
	}


    void resize(int _screen_width, int _screen_height) {
        // make sure the viewport matches the new window dimensions; note that width and 
        // height will be significantly larger than specified on retina displays.
        glViewport(0, 0, _screen_width, _screen_height);
        screen_width = _screen_width;
        screen_height = _screen_height;
        updateUboData();
        
    }

    // render loop
	void render(const Scene & _scene, bool render_ImGUI = true) {
        // update the UBO data
        updateUboData();

        // upload the UBO data
        uploadUbo();
        
        
        
        // clear screen
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        

        // get light info
        const std::vector<std::unique_ptr<Light>> & sceneLights = _scene.sceneLights;
        // currently we only support one light and it is a point light
        if (sceneLights.size() > 0) {
            Light * light = sceneLights[0].get();
            if (light->lightType == POINT_LIGHT) {
                PointLight * pointLight = dynamic_cast<PointLight*>(light);
                // do something with the point light, but not implemented yet...
            }
        }

        // compute pass
        const std::vector<std::shared_ptr<ComputeComponent>> & computeQueue = _scene.computeQueue;
        // compute the scene using openGL compute shader
        for (auto computeComponent : computeQueue) {
            if (computeComponent->active) {
                computeComponent->Compute();
            }
        }

        // TODO: probably need to implement a render pass logic here
        // so that user can define the render pass order and the framebuffer will be automatically generated and switched
        // glBindFramebuffer(GL_FRAMEBUFFER, 1); // use the default framebuffer
        // render pass
        const std::vector<std::shared_ptr<RenderComponent>> & renderQueue = _scene.renderQueue;
        // render the scene using openGL rasterization pipeline
        for (auto renderComponent : renderQueue) {
            if (renderComponent->active) {
                renderComponent->Render();
            }
        }
        

	}


    private:
    //----------------------------------------------------------------------------------
    const bool _vSync = true; // Enable vsync

};












#endif