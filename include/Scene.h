#ifndef SCENE_H
#define SCENE_H




#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>

#include "CommonSceneObject.h"
#include <Shader.h>

#include <iostream>
#include <chrono>

#include "DestructiveCSSceneObject.h"

enum SceneName {
    SCENE1 = 0,
    SCENE2 = 1,
    SCENE3 = 2,
    SCENE4 = 3,
};


class Scene {
public:
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
    }


    Scene(RenderContext* _renderContext = nullptr, SceneName sceneName = SCENE1) {
        switch (sceneName) {
            case SCENE1:
                createScene1(_renderContext);
                break;
            case SCENE2:
                createScene2(_renderContext);
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
        GModel * armadillo = new GModel("resource/stanford-bunny1x1x1.obj");
        CommonSceneObject * CommonSceneObject_ARMADILLO = new CommonSceneObject(armadillo,TRANSPARENT,renderContext);
        CommonSceneObject_ARMADILLO->setMaterial(LAMBERTIAN, 0.0, glm::vec4(0.6, 0.7, 1.0, 0.4));
        CommonSceneObject_ARMADILLO->setModelMatrix(glm::translate(glm::mat4(1.0), glm::vec3(0, 0, 0)) * glm::scale(glm::mat4(1.0), glm::vec3(0.91, 0.91, 0.91)));
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

    // old scene from my previous project RTRT， just for reference
    bool createScene2(RenderContext* _renderContext){
//        renderContext = _renderContext;
//    
//        glm::vec3 rotationAxis = glm::vec3(0, 1, 0);
//        glm::mat4 identity = glm::mat4(1.0);
//        glm::vec3 v0(-1, 2, -0.2);
//        glm::vec3 v1(1, 2, 0.2);
//        glm::vec3 v2(0, 4, 0);
//        
//        
//        // skybox
//        skyboxTexture = new SkyboxTexture();
//        skyboxTexture->loadFromFolder("resource/skybox");
//        skyboxTexture->createSkyboxTexture();
//        GSkybox * _skybox = new GSkybox();
//        _skybox->setShader(new Shader("shaders/skybox_shader.vert", "shaders/skybox_shader.frag"));
//        _skybox->setTexture(skyboxTexture);
//        Skybox * skybox = new Skybox(_skybox);
//        skybox->attachToSceneRenderList(renderQueue);
//        sceneObjects.push_back(skybox);
//        
//        // waifu, would be a performance bottleneck because of the high triangle count
//        Texture * waifuTexture = new Texture();
//        waifuTexture->loadFromFile("resource/mebius_diffuse.png");
//        waifuTexture->removeAlphaChannel();
//        waifuTexture->resizeData(1024, 1024,3);
//        waifuTexture->createGPUTexture();
//        GModel* waifu = new GModel("resource/mebius.obj");
//        waifu->setSkyboxTexture(skyboxTexture);
//        CommonSceneObject * CommonSceneObject7 = new CommonSceneObject(waifu,OPAQUE,renderContext);
//        CommonSceneObject7->setMaterial(LAMBERTIAN, 1.5, glm::vec4(1.0,1.0,1.0, 1.0), waifuTexture);
//        CommonSceneObject7->setModelMatrix(glm::translate(glm::mat4(1.0), glm::vec3(0, 2, 2.2)) * glm::rotate(identity, glm::radians(-95.0f), rotationAxis) * glm::scale(glm::mat4(1.0), glm::vec3(0.08, 0.08, 0.08)));
//        CommonSceneObject7->attachToSceneRenderList(renderQueue);
//        sceneObjects.push_back(CommonSceneObject7);
//        
//        
//        // the ground sphere
//        GSphere * sphereObject2 = new GSphere();
//        sphereObject2->setSkyboxTexture(skyboxTexture);
//        CommonSceneObject * CommonSceneObject2 = new CommonSceneObject(sphereObject2,OPAQUE,renderContext);
//        CommonSceneObject2->setMaterial(LAMBERTIAN, 0.0, glm::vec4(0.5, 0.5, 0.5, 1.0));
//        CommonSceneObject2->setModelMatrix(glm::translate(glm::mat4(1.0), glm::vec3(0, -100, 0)) * glm::scale(glm::mat4(1.0), glm::vec3(100, 100, 100)));
//        CommonSceneObject2->attachToSceneRenderList(renderQueue);
//        sceneObjects.push_back(CommonSceneObject2);
//        
//        // the earth sphere
//        Texture * earthTexture = new Texture();
//        earthTexture->loadFromFile("resource/earthmap.jpg");
//        earthTexture->createGPUTexture();
//        GSphere* sphereObject3 = new GSphere();
//        sphereObject3->setSkyboxTexture(skyboxTexture);
//        CommonSceneObject * CommonSceneObject3 = new CommonSceneObject(sphereObject3,OPAQUE,renderContext);
//        CommonSceneObject3->setMaterial(LAMBERTIAN, 0.0, glm::vec4(1.0, 1.0, 1.0, 1.0), earthTexture);
//        CommonSceneObject3->setModelMatrix(glm::translate(glm::mat4(1.0), glm::vec3(0, 1, 2.2)) * glm::rotate(identity, glm::radians(90.0f), rotationAxis) *glm::scale(glm::mat4(1.0), glm::vec3(1, 1, 1)));
//        CommonSceneObject3->attachToSceneRenderList(renderQueue);
//        sceneObjects.push_back(CommonSceneObject3);
//        
//        // the metal sphere
//        GSphere* sphereObject4 = new GSphere();
//        sphereObject4->setSkyboxTexture(skyboxTexture);
//        CommonSceneObject * CommonSceneObject4 = new CommonSceneObject(sphereObject4,OPAQUE,renderContext);
//        CommonSceneObject4->setMaterial(METAL, 0.2, glm::vec4(0.7, 0.6, 0.5, 1.0));
//        CommonSceneObject4->setModelMatrix(glm::translate(glm::mat4(1.0), glm::vec3(0, 1, 0)) * glm::scale(glm::mat4(1.0), glm::vec3(1, 1, 1)));
//        CommonSceneObject4->attachToSceneRenderList(renderQueue);
//        sceneObjects.push_back(CommonSceneObject4);
//
//
//        // diffuse light sphere
//        GSphere* sphereObject5 = new GSphere();
//        sphereObject5->setSkyboxTexture(skyboxTexture);
//        CommonSceneObject * CommonSceneObject5 = new CommonSceneObject(sphereObject5,OPAQUE,renderContext);
//        CommonSceneObject5->setMaterial(EMISSIVE, 0.0, glm::vec4(2, 2, 2, 1.0));
//        CommonSceneObject5->setModelMatrix(glm::translate(glm::mat4(1.0), glm::vec3(1.5, 0.45, 0)) * glm::scale(glm::mat4(1.0), glm::vec3(0.5, 0.5, 0.5)));
//        CommonSceneObject5->attachToSceneRenderList(renderQueue);
//        sceneObjects.push_back(CommonSceneObject5);
//
//        // record the light source to sceneObjects
//        PointLight * sphereObject5Light = new PointLight(glm::vec3(1.5, 0.45, 0), glm::vec3(2, 2, 2));
//        sceneObjects.push_back(sphereObject5Light); // then such light can be updated in the logic loop by ticking its component
//        sceneLights.push_back(sphereObject5Light);
//        // add periodic translation component to the light source mesh
//        std::unique_ptr<ObjectPeriodicTranslationComponent> sphereObject5TranslationComponent = std::make_unique<ObjectPeriodicTranslationComponent>(glm::vec3(0, 0, 1), glm::vec3(0, 0, 1));
//        CommonSceneObject5->addComponent(std::move(sphereObject5TranslationComponent));
//        // we also update the light source Light object in the logic loop by ticking its component
//        std::unique_ptr<ObjectPeriodicTranslationComponent> sphereObject5TranslationComponent2 = std::make_unique<ObjectPeriodicTranslationComponent>(glm::vec3(0, 0, 1), glm::vec3(0, 0, 1));
//        sphereObject5Light->addComponent(std::move(sphereObject5TranslationComponent2));
//        
//        
//        GTriangle* triangleObject = new GTriangle(v0, v1, v2);
//        triangleObject->setSkyboxTexture(skyboxTexture);
//        CommonSceneObject * CommonSceneObject6 = new CommonSceneObject(triangleObject,OPAQUE,renderContext);
//        CommonSceneObject6->setMaterial(LAMBERTIAN, 0.0, glm::vec4(0.5, 0.5, 0.5, 1.0));
//        CommonSceneObject6->setModelMatrix(glm::translate(glm::mat4(1.0), glm::vec3(0, 0, 0)) * glm::scale(glm::mat4(1.0), glm::vec3(1, 1, 1)));
//        CommonSceneObject6->attachToSceneRenderList(renderQueue);
//        sceneObjects.push_back(CommonSceneObject6);
//        
//        // the cube
//        Texture * cubeTexture = new Texture();
//        cubeTexture->loadFromFile("resource/night.png");
//        cubeTexture->createGPUTexture();
//        GModel* cubeObject = new GModel("resource/cube.obj");
//        cubeObject->setSkyboxTexture(skyboxTexture);
//        CommonSceneObject * CommonSceneObject8 = new CommonSceneObject(cubeObject,OPAQUE,renderContext);
//        CommonSceneObject8->setMaterial(LAMBERTIAN, 0.0, glm::vec4(1.0, 1.0, 1.0, 1.0), cubeTexture);
//        CommonSceneObject8->setModelMatrix(glm::translate(glm::mat4(1.0), glm::vec3(1.5, 0.5, 2.0)) * glm::scale(glm::mat4(1.0), glm::vec3(0.5, 0.5, 0.5)));
//        CommonSceneObject8->attachToSceneRenderList(renderQueue);
//        // add rotation component to the cube
//        std::unique_ptr<ObjectRotationComponent> cubeRotationComponent = std::make_unique<ObjectRotationComponent>(4.0f);
//        CommonSceneObject8->addComponent(std::move(cubeRotationComponent));
//        sceneObjects.push_back(CommonSceneObject8);
//
//        
//        // the glass sphere (transparent object in rendering order)
//        GSphere * sphereObject = new GSphere();
//        sphereObject->setSkyboxTexture(skyboxTexture);
//        CommonSceneObject * CommonSceneObject1 = new CommonSceneObject(sphereObject,TRANSPARENT,renderContext);
//        CommonSceneObject1->setMaterial(DIELECTRIC, 1.5, glm::vec4(0.3, 0.4, 0.8, 0.6));
//        CommonSceneObject1->setModelMatrix(glm::translate(glm::mat4(1.0), glm::vec3(0, 1, -2.2)) * glm::scale(glm::mat4(1.0), glm::vec3(1, 1, 1)));
//        CommonSceneObject1->attachToSceneRenderList(renderQueue);
//        sceneObjects.push_back(CommonSceneObject1);
//        return true;
        return false;
    }

};






#endif