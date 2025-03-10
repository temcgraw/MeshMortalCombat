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
    // TODO: use smart pointers
    // and fix all memory leaks...
    std::vector<SceneObject*> sceneObjects;
    std::vector<std::shared_ptr<RenderComponent>> renderQueue; // object components to be rendered each frame
    std::vector<std::shared_ptr<ComputeComponent>> computeQueue; // object components to be computed via compute shader each frame
    std::vector<Light*> sceneLights; // point lights in the scene
    
    SkyboxTexture * skyboxTexture; // skybox texture
    RenderContext * renderContext = nullptr;

    void sortRenderQueue() {
        std::sort(renderQueue.begin(), renderQueue.end(), [](std::shared_ptr<RenderComponent> a, std::shared_ptr<RenderComponent> b) {
            return a->renderPriority < b->renderPriority;
        });
    }

    Scene() {
    }

    virtual ~Scene() {
        for (SceneObject* sceneObject : sceneObjects) {
            delete sceneObject;
        }
        for (Light* light : sceneLights) {
            delete light;
        }
        
    }


    Scene(RenderContext* _renderContext = nullptr, SceneName sceneName = SCENE1) {
        switch (sceneName) {
            case SCENE1:
                createScene1(_renderContext);
                break;
            case SCENE2:
                createScene2(_renderContext);
                break;
            case SCENE3:
                createScene3(_renderContext);
                break;
            case SCENE4:
                createScene4(_renderContext);
                break;
            case SCENE5:
                createScene5(_renderContext);
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

    // the Destructible Adaptive Grid Particle System scene
    bool createScene1(RenderContext* _renderContext){
        renderContext = _renderContext;
        // skybox
        skyboxTexture = new SkyboxTexture();
        skyboxTexture->loadFromFolder("resource/skybox");
        skyboxTexture->createSkyboxTexture();
        GSkybox * _skybox = new GSkybox();
        _skybox->setShader(new Shader("shaders/skybox_shader.vert", "shaders/skybox_shader.frag"));
        _skybox->setTexture(skyboxTexture);
        Skybox * skybox = new Skybox(_skybox);
        skybox->attachToSceneRenderList(renderQueue);
        sceneObjects.push_back(skybox);
        // --------------------------------------------------------------------------------
        // probably need to refactor this part, it's not a good idea to put one object into another object
        // maybe use GModel instead of CommonSceneObject
        // the armadillo we want to voxelize and destruct
        GModel * armadillo = new GModel("resource/armtest1x1x1.obj");
        CommonSceneObject * CommonSceneObject_ARMADILLO = new CommonSceneObject(armadillo,TRANSPARENT,renderContext);
        CommonSceneObject_ARMADILLO->setMaterial(LAMBERTIAN, 0.0, glm::vec4(0.6, 0.7, 1.0, 0.4));
        CommonSceneObject_ARMADILLO->setModelMatrix(glm::translate(glm::mat4(1.0), glm::vec3(-0.5, -0.5, -0.5)) * glm::scale(glm::mat4(1.0), glm::vec3(1,1,1)));
        // we don't want it to individually render, so we don't attach it to the renderQueue
        // but instead, it's rendering is controlled by the DestructiveCSSceneObject
        //CommonSceneObject_ARMADILLO->attachToSceneRenderList(renderQueue);
        sceneObjects.push_back(CommonSceneObject_ARMADILLO);
        
        // the DestructiveCSSceneObject
        DestructiveCSSceneObject * destructiveCSSceneObject = new DestructiveCSSceneObject(std::shared_ptr<GModel>(armadillo), TRANSPARENT, renderContext);
        destructiveCSSceneObject->attachToSceneComputeList(computeQueue);
        destructiveCSSceneObject->attachToSceneRenderList(renderQueue);
        sceneObjects.push_back(destructiveCSSceneObject);


        return true;
    }
    bool createScene2(RenderContext* _renderContext){
        renderContext = _renderContext;
        skyboxTexture = new SkyboxTexture();
        skyboxTexture->loadFromFolder("resource/skybox");
        skyboxTexture->createSkyboxTexture();
        GSkybox * _skybox = new GSkybox();
        _skybox->setShader(new Shader("shaders/skybox_shader.vert", "shaders/skybox_shader.frag"));
        _skybox->setTexture(skyboxTexture);
        Skybox * skybox = new Skybox(_skybox);
        skybox->attachToSceneRenderList(renderQueue);
        sceneObjects.push_back(skybox);
        GModel * armadillo = new GModel("resource/stanford-bunny1x1x1.obj");
        CommonSceneObject * CommonSceneObject_ARMADILLO = new CommonSceneObject(armadillo,TRANSPARENT,renderContext);
        CommonSceneObject_ARMADILLO->setMaterial(LAMBERTIAN, 0.0, glm::vec4(0.6, 0.7, 1.0, 0.4));
        CommonSceneObject_ARMADILLO->setModelMatrix(glm::translate(glm::mat4(1.0), glm::vec3(-0.5, -0.5, -0.5)) * glm::scale(glm::mat4(1.0), glm::vec3(1,1,1)));
        sceneObjects.push_back(CommonSceneObject_ARMADILLO);
        DestructiveCSSceneObject * destructiveCSSceneObject = new DestructiveCSSceneObject(std::shared_ptr<GModel>(armadillo), TRANSPARENT, renderContext);
        destructiveCSSceneObject->attachToSceneComputeList(computeQueue);
        destructiveCSSceneObject->attachToSceneRenderList(renderQueue);
        sceneObjects.push_back(destructiveCSSceneObject);
        return true;
    }
    bool createScene3(RenderContext* _renderContext){
        renderContext = _renderContext;
        skyboxTexture = new SkyboxTexture();
        skyboxTexture->loadFromFolder("resource/skybox");
        skyboxTexture->createSkyboxTexture();
        GSkybox * _skybox = new GSkybox();
        _skybox->setShader(new Shader("shaders/skybox_shader.vert", "shaders/skybox_shader.frag"));
        _skybox->setTexture(skyboxTexture);
        Skybox * skybox = new Skybox(_skybox);
        skybox->attachToSceneRenderList(renderQueue);
        sceneObjects.push_back(skybox);
        GModel * armadillo = new GModel("resource/gummyBear1x1x1.obj");
        CommonSceneObject * CommonSceneObject_ARMADILLO = new CommonSceneObject(armadillo,TRANSPARENT,renderContext);
        CommonSceneObject_ARMADILLO->setMaterial(LAMBERTIAN, 0.0, glm::vec4(0.6, 0.7, 1.0, 0.4));
        CommonSceneObject_ARMADILLO->setModelMatrix(glm::translate(glm::mat4(1.0), glm::vec3(-0.5, -0.5, -0.5)) * glm::scale(glm::mat4(1.0), glm::vec3(1,1,1)));
        sceneObjects.push_back(CommonSceneObject_ARMADILLO);
        DestructiveCSSceneObject * destructiveCSSceneObject = new DestructiveCSSceneObject(std::shared_ptr<GModel>(armadillo), TRANSPARENT, renderContext);
        destructiveCSSceneObject->attachToSceneComputeList(computeQueue);
        destructiveCSSceneObject->attachToSceneRenderList(renderQueue);
        sceneObjects.push_back(destructiveCSSceneObject);
        return true;
    }
    bool createScene4(RenderContext* _renderContext){
        renderContext = _renderContext;
        skyboxTexture = new SkyboxTexture();
        skyboxTexture->loadFromFolder("resource/skybox");
        skyboxTexture->createSkyboxTexture();
        GSkybox * _skybox = new GSkybox();
        _skybox->setShader(new Shader("shaders/skybox_shader.vert", "shaders/skybox_shader.frag"));
        _skybox->setTexture(skyboxTexture);
        Skybox * skybox = new Skybox(_skybox);
        skybox->attachToSceneRenderList(renderQueue);
        sceneObjects.push_back(skybox);
        GModel * armadillo = new GModel("resource/skyscraper1x1x1.obj");
        CommonSceneObject * CommonSceneObject_ARMADILLO = new CommonSceneObject(armadillo,TRANSPARENT,renderContext);
        CommonSceneObject_ARMADILLO->setMaterial(LAMBERTIAN, 0.0, glm::vec4(0.6, 0.7, 1.0, 0.4));
        CommonSceneObject_ARMADILLO->setModelMatrix(glm::translate(glm::mat4(1.0), glm::vec3(-0.5, -0.5, -0.5)) * glm::scale(glm::mat4(1.0), glm::vec3(1,1,1)));
        sceneObjects.push_back(CommonSceneObject_ARMADILLO);
        DestructiveCSSceneObject * destructiveCSSceneObject = new DestructiveCSSceneObject(std::shared_ptr<GModel>(armadillo), TRANSPARENT, renderContext);
        destructiveCSSceneObject->attachToSceneComputeList(computeQueue);
        destructiveCSSceneObject->attachToSceneRenderList(renderQueue);
        sceneObjects.push_back(destructiveCSSceneObject);
        return true;
    }
    bool createScene5(RenderContext* _renderContext){
        renderContext = _renderContext;
        skyboxTexture = new SkyboxTexture();
        skyboxTexture->loadFromFolder("resource/skybox");
        skyboxTexture->createSkyboxTexture();
        GSkybox * _skybox = new GSkybox();
        _skybox->setShader(new Shader("shaders/skybox_shader.vert", "shaders/skybox_shader.frag"));
        _skybox->setTexture(skyboxTexture);
        Skybox * skybox = new Skybox(_skybox);
        skybox->attachToSceneRenderList(renderQueue);
        sceneObjects.push_back(skybox);
        GModel * armadillo = new GModel("resource/Teapot1x1x1.obj");
        CommonSceneObject * CommonSceneObject_ARMADILLO = new CommonSceneObject(armadillo,TRANSPARENT,renderContext);
        CommonSceneObject_ARMADILLO->setMaterial(LAMBERTIAN, 0.0, glm::vec4(0.6, 0.7, 1.0, 0.4));
        CommonSceneObject_ARMADILLO->setModelMatrix(glm::translate(glm::mat4(1.0), glm::vec3(-0.5, -0.5, -0.5)) * glm::scale(glm::mat4(1.0), glm::vec3(1,1,1)));
        sceneObjects.push_back(CommonSceneObject_ARMADILLO);
        DestructiveCSSceneObject * destructiveCSSceneObject = new DestructiveCSSceneObject(std::shared_ptr<GModel>(armadillo), TRANSPARENT, renderContext);
        destructiveCSSceneObject->attachToSceneComputeList(computeQueue);
        destructiveCSSceneObject->attachToSceneRenderList(renderQueue);
        sceneObjects.push_back(destructiveCSSceneObject);
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
        scene = new Scene(renderContext, sceneName);
    }
    void setRenderContext(RenderContext* _renderContext) {
        renderContext = _renderContext;
    }
    Scene* getScene() {
        return scene;
    }
    
    std::vector<SceneObject*> & getSceneObjects() {
        if(scene == nullptr) {
            std::cout<<"[SceneManager]: scene is nullptr"<<std::endl;
            static std::vector<SceneObject*> emptyVector;
            return emptyVector;
        }
        return scene->sceneObjects;
    }

private:
    Scene* scene;// TODO; use smart pointer
    RenderContext* renderContext = nullptr;
    SceneManager(){}
    SceneManager(const SceneManager&) = delete;
    SceneManager& operator=(const SceneManager&) = delete;

};


#endif