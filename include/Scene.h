#ifndef SCENE_H
#define SCENE_H




#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>

#include "CommonSceneObject.h"
#include "DestructiveCSSceneObject.h"// TODO: decouple this, use other method to create scene with DestructiveCSSceneObject
#include <Shader.h>

#include <iostream>
#include <chrono>



enum SceneName {
    SCENE1 = 0,
    SCENE2 = 1,
    SCENE3 = 2,
    SCENE4 = 3,
    SCENE5 = 4
};


class Scene {
public:
    std::vector<std::unique_ptr<SceneObject>> sceneObjects;
    std::vector<std::shared_ptr<RenderComponent>> renderQueue; // object components to be rendered each frame
    std::vector<std::shared_ptr<ComputeComponent>> computeQueue; // object components to be computed via compute shader each frame
    std::vector<std::unique_ptr<Light>> sceneLights;
    std::shared_ptr<SkyboxTexture> skyboxTexture;
    std::vector<std::unique_ptr<Shader>> shaders; // scene will manage the shaders, when scene is destroyed, the shaders will be destroyed

    void sortRenderQueue() {
        std::sort(renderQueue.begin(), renderQueue.end(), [](std::shared_ptr<RenderComponent> a, std::shared_ptr<RenderComponent> b) {
            return a->renderPriority < b->renderPriority;
        });
    }

    Scene() {
    }

    virtual ~Scene() {
    }


    Scene(SceneName sceneName = SCENE1) {
        switch (sceneName) {
            case SCENE1:
                createScene1();
                break;
            case SCENE2:
                createScene2();
                break;
            case SCENE3:
                createScene3();
                break;
            case SCENE4:
                createScene4();
                break;
            case SCENE5:
                createScene5();
                break;
            default:
                std::cout<<"[Scene]: invalid scene name"<<std::endl;
                break;
        }
        
        
        // sort the renderQueue based on renderPriority
        // the lower the renderPriority, the earlier it is rendered
        std::sort(renderQueue.begin(), renderQueue.end(), [](std::shared_ptr<RenderComponent> a, std::shared_ptr<RenderComponent> b) {
            return a->renderPriority < b->renderPriority;
        });



    
    }

    const char * vertexShaderPath = "shaders/object_shader.vert";
    const char * fragmentShaderPath = "shaders/object_shader.frag";

    // the Destructible Adaptive Grid Particle System scene
    bool createScene1(){
        // create skybox texture
        skyboxTexture = std::make_shared<SkyboxTexture>();
        skyboxTexture->loadFromFolder("resource/skybox");
        skyboxTexture->createSkyboxTexture();
        // shaders are managed by the scene, so we use unique_ptr
        auto skyboxShader = std::make_unique<Shader>("shaders/skybox_shader.vert", "shaders/skybox_shader.frag");
        // get the raw pointer of the shader, and pass it to the skybox object
        Shader* pSkyboxShader = skyboxShader.get();
        shaders.push_back(std::move(skyboxShader));
        // create the skybox GObject and SceneObject
        auto _skybox = std::make_unique<GSkybox>();
        _skybox->setShader(pSkyboxShader);
        _skybox->setTexture(skyboxTexture.get());
        auto skybox = std::make_unique<Skybox>(std::move(_skybox));
        skybox->attachToSceneRenderList(renderQueue);
        sceneObjects.push_back(std::move(skybox));
        auto objectShader = std::make_unique<Shader>("shaders/object_shader.vert", "shaders/object_shader.frag");
        Shader* pObjectShader = objectShader.get();
        shaders.push_back(std::move(objectShader));
        auto inputMesh = std::make_unique<GModel>("resource/stanford-bunny1x1x1.obj");
        inputMesh->setModel(glm::scale(glm::mat4(1.0f), glm::vec3(1.25f)) * glm::translate(glm::mat4(1.0f), glm::vec3(-0.5f, -0.25f, -0.5f)));
        inputMesh->setShader(pObjectShader);
        inputMesh->setColor(glm::vec4(0.6, 0.7, 1.0, 0.4));
        auto destructiveCSSceneObject = std::make_unique<DestructiveCSSceneObject>(std::shared_ptr<GModel>(std::move(inputMesh)), TRANSPARENT);
        destructiveCSSceneObject->attachToSceneComputeList(computeQueue);
        destructiveCSSceneObject->attachToSceneRenderList(renderQueue);
        sceneObjects.push_back(std::move(destructiveCSSceneObject));

        return true;
    }
    bool createScene2(){
        // create skybox texture
        skyboxTexture = std::make_shared<SkyboxTexture>();
        skyboxTexture->loadFromFolder("resource/skybox");
        skyboxTexture->createSkyboxTexture();
        // shaders are managed by the scene, so we use unique_ptr
        auto skyboxShader = std::make_unique<Shader>("shaders/skybox_shader.vert", "shaders/skybox_shader.frag");
        // get the raw pointer of the shader, and pass it to the skybox object
        Shader* pSkyboxShader = skyboxShader.get();
        shaders.push_back(std::move(skyboxShader));
        // create the skybox GObject and SceneObject
        auto _skybox = std::make_unique<GSkybox>();
        _skybox->setShader(pSkyboxShader);
        _skybox->setTexture(skyboxTexture.get());
        auto skybox = std::make_unique<Skybox>(std::move(_skybox));
        skybox->attachToSceneRenderList(renderQueue);
        sceneObjects.push_back(std::move(skybox));
        auto objectShader = std::make_unique<Shader>("shaders/object_shader.vert", "shaders/object_shader.frag");
        Shader* pObjectShader = objectShader.get();
        shaders.push_back(std::move(objectShader));
        auto inputMesh = std::make_unique<GModel>("resource/armtest1x1x1.obj");
        inputMesh->setModel(glm::scale(glm::mat4(1.0f), glm::vec3(1.25f)) * glm::translate(glm::mat4(1.0f), glm::vec3(-0.5f, -0.25f, -0.5f)));
        inputMesh->setShader(pObjectShader);
        inputMesh->setColor(glm::vec4(0.6, 0.7, 1.0, 0.4));
        auto destructiveCSSceneObject = std::make_unique<DestructiveCSSceneObject>(std::shared_ptr<GModel>(std::move(inputMesh)), TRANSPARENT);
        destructiveCSSceneObject->attachToSceneComputeList(computeQueue);
        destructiveCSSceneObject->attachToSceneRenderList(renderQueue);
        sceneObjects.push_back(std::move(destructiveCSSceneObject));

        return true;
    }
    bool createScene3(){
        // create skybox texture
        skyboxTexture = std::make_shared<SkyboxTexture>();
        skyboxTexture->loadFromFolder("resource/skybox");
        skyboxTexture->createSkyboxTexture();
        // shaders are managed by the scene, so we use unique_ptr
        auto skyboxShader = std::make_unique<Shader>("shaders/skybox_shader.vert", "shaders/skybox_shader.frag");
        // get the raw pointer of the shader, and pass it to the skybox object
        Shader* pSkyboxShader = skyboxShader.get();
        shaders.push_back(std::move(skyboxShader));
        // create the skybox GObject and SceneObject
        auto _skybox = std::make_unique<GSkybox>();
        _skybox->setShader(pSkyboxShader);
        _skybox->setTexture(skyboxTexture.get());
        auto skybox = std::make_unique<Skybox>(std::move(_skybox));
        skybox->attachToSceneRenderList(renderQueue);
        sceneObjects.push_back(std::move(skybox));
        auto objectShader = std::make_unique<Shader>("shaders/object_shader.vert", "shaders/object_shader.frag");
        Shader* pObjectShader = objectShader.get();
        shaders.push_back(std::move(objectShader));
        auto inputMesh = std::make_unique<GModel>("resource/gummyBear1x1x1.obj");
        inputMesh->setModel(glm::scale(glm::mat4(1.0f), glm::vec3(1.25f)) * glm::translate(glm::mat4(1.0f), glm::vec3(-0.5f, -0.25f, -0.5f)));
        inputMesh->setShader(pObjectShader);
        inputMesh->setColor(glm::vec4(0.6, 0.7, 1.0, 0.4));
        auto destructiveCSSceneObject = std::make_unique<DestructiveCSSceneObject>(std::shared_ptr<GModel>(std::move(inputMesh)), TRANSPARENT);
        destructiveCSSceneObject->attachToSceneComputeList(computeQueue);
        destructiveCSSceneObject->attachToSceneRenderList(renderQueue);
        sceneObjects.push_back(std::move(destructiveCSSceneObject));

        return true;
    }
    bool createScene4(){
        // create skybox texture
        skyboxTexture = std::make_shared<SkyboxTexture>();
        skyboxTexture->loadFromFolder("resource/skybox");
        skyboxTexture->createSkyboxTexture();
        // shaders are managed by the scene, so we use unique_ptr
        auto skyboxShader = std::make_unique<Shader>("shaders/skybox_shader.vert", "shaders/skybox_shader.frag");
        // get the raw pointer of the shader, and pass it to the skybox object
        Shader* pSkyboxShader = skyboxShader.get();
        shaders.push_back(std::move(skyboxShader));
        // create the skybox GObject and SceneObject
        auto _skybox = std::make_unique<GSkybox>();
        _skybox->setShader(pSkyboxShader);
        _skybox->setTexture(skyboxTexture.get());
        auto skybox = std::make_unique<Skybox>(std::move(_skybox));
        skybox->attachToSceneRenderList(renderQueue);
        sceneObjects.push_back(std::move(skybox));
        auto objectShader = std::make_unique<Shader>("shaders/object_shader.vert", "shaders/object_shader.frag");
        Shader* pObjectShader = objectShader.get();
        shaders.push_back(std::move(objectShader));
        auto inputMesh = std::make_unique<GModel>("resource/skyscraper1x1x1.obj");
        inputMesh->setModel(glm::scale(glm::mat4(1.0f), glm::vec3(1.25f)) * glm::translate(glm::mat4(1.0f), glm::vec3(-0.5f, -0.25f, -0.5f)));
        inputMesh->setShader(pObjectShader);
        inputMesh->setColor(glm::vec4(0.6, 0.7, 1.0, 0.4));
        auto destructiveCSSceneObject = std::make_unique<DestructiveCSSceneObject>(std::shared_ptr<GModel>(std::move(inputMesh)), TRANSPARENT);
        destructiveCSSceneObject->attachToSceneComputeList(computeQueue);
        destructiveCSSceneObject->attachToSceneRenderList(renderQueue);
        sceneObjects.push_back(std::move(destructiveCSSceneObject));

        return true;
    }
    bool createScene5(){
        // create skybox texture
        skyboxTexture = std::make_shared<SkyboxTexture>();
        skyboxTexture->loadFromFolder("resource/skybox");
        skyboxTexture->createSkyboxTexture();
        // shaders are managed by the scene, so we use unique_ptr
        auto skyboxShader = std::make_unique<Shader>("shaders/skybox_shader.vert", "shaders/skybox_shader.frag");
        // get the raw pointer of the shader, and pass it to the skybox object
        Shader* pSkyboxShader = skyboxShader.get();
        shaders.push_back(std::move(skyboxShader));
        // create the skybox GObject and SceneObject
        auto _skybox = std::make_unique<GSkybox>();
        _skybox->setShader(pSkyboxShader);
        _skybox->setTexture(skyboxTexture.get());
        auto skybox = std::make_unique<Skybox>(std::move(_skybox));
        skybox->attachToSceneRenderList(renderQueue);
        sceneObjects.push_back(std::move(skybox));
        auto objectShader = std::make_unique<Shader>("shaders/object_shader.vert", "shaders/object_shader.frag");
        Shader* pObjectShader = objectShader.get();
        shaders.push_back(std::move(objectShader));
        auto inputMesh = std::make_unique<GModel>("resource/Teapot1x1x1.obj");
        inputMesh->setModel(glm::scale(glm::mat4(1.0f), glm::vec3(1.25f)) * glm::translate(glm::mat4(1.0f), glm::vec3(-0.5f, -0.25f, -0.5f)));
        inputMesh->setShader(pObjectShader);
        inputMesh->setColor(glm::vec4(0.6, 0.7, 1.0, 0.4));
        auto destructiveCSSceneObject = std::make_unique<DestructiveCSSceneObject>(std::shared_ptr<GModel>(std::move(inputMesh)), TRANSPARENT);
        destructiveCSSceneObject->attachToSceneComputeList(computeQueue);
        destructiveCSSceneObject->attachToSceneRenderList(renderQueue);
        sceneObjects.push_back(std::move(destructiveCSSceneObject));

        return true;
    }

};

// the application can only have one scene manager and one scene active at a time
// so we use singleton pattern here
class SceneManager {
public:
    static SceneManager& GetInstance() {
        static SceneManager instance;
        return instance;
    }
    void switchScene(SceneName sceneName = SCENE1) {
        if (scene) {
            delete scene;
        }
        scene = new Scene(sceneName);
    }

    Scene* getScene() {
        return scene;
    }
    
    std::vector<std::unique_ptr<SceneObject>> & getSceneObjects() {
        if(scene == nullptr) {
            std::cout<<"[SceneManager]: scene is nullptr"<<std::endl;
            static std::vector<std::unique_ptr<SceneObject>> emptyVector;
            return emptyVector;
        }
        return scene->sceneObjects;
    }

private:
    Scene* scene;// TODO; use smart pointer
    SceneManager(){}
    SceneManager(const SceneManager&) = delete;
    SceneManager& operator=(const SceneManager&) = delete;

};


#endif