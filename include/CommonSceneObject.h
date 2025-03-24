#ifndef CommonSceneObject_H
#define CommonSceneObject_H


#include "GraphicObject.h"
#include "SceneObject.h"

#include <memory>

enum MaterialType {
    LAMBERTIAN = 0,
    METAL = 1,
    DIELECTRIC = 2,
    EMISSIVE = 3
};

class ObjectRenderComponent : public RenderComponent {
public:
    ObjectRenderComponent(std::shared_ptr<GMVPObject> _obj, enum renderQueue _renderPriority = OPAQUE) {
        this->obj = _obj;
        renderPriority = _renderPriority;
    }
    void Render() override {

        obj->draw();
    }

private:
std::shared_ptr<GMVPObject> obj;
};





// CommonSceneObject class, it is designed to combine the GObject and RenderableSceneObject together
// GObject is for low-level rendering data and pipeline management, while RenderableSceneObject is for high-level scene management
// instead of inheriting from GObject class, it contains an GObject pointer, and it can be constructed from an GObject pointer
// the key data of this class is the model matrix, material type, and texture, making it easy to modify the rendering data and
// interact with other components in CPU
class CommonSceneObject : public RenderableSceneObject {
public:
    CommonSceneObject(std::shared_ptr<GMVPObject> _obj, enum renderQueue _renderPriority = OPAQUE) {
        this->obj = _obj;
        this->renderComponent = std::make_shared<ObjectRenderComponent>(obj,_renderPriority);
        
    }
    ~CommonSceneObject() {
    }

    void setModelMatrix(glm::mat4 modelMatrix) {
        obj->setModel(modelMatrix);
        this->modelMatrix = modelMatrix;
    }

    glm::mat4 getModelMatrix() {
        return modelMatrix;
    }

    void setMaterial(int type, float fuzzOrIOR, glm::vec4 baseColor, Texture * _texture = nullptr) {
        // ensure if fuzziness, it is in the range [0, 1]
        materialType = type;
        obj->setColor(baseColor);
        // we specify the corresponding shader in default renderer to make it easy to manage the mapping between raytracing material type and rasterization shading shader 
        const char * vertexShaderPath = "shaders/object_shader.vert";
        const char * fragmentShaderPath = "shaders/object_shader.frag";
        if (_texture != nullptr) {
            obj->setTexture(_texture);
            this->texture = _texture;
            // resize the texture to 1024x1024x3
            texture->resizeTexture(1024, 1024,3);
        }
        if (type == 0) {
            obj->setShader(new Shader(vertexShaderPath, fragmentShaderPath));
        }
        else if (type == 1) {
            obj->setShader(new Shader(vertexShaderPath, fragmentShaderPath));
            // hard-coded metal material uniforms
            obj->shader->use();
            obj->shader->setFloat("averageSlope", std::min(std::max(fuzzOrIOR, 0.0f), 1.0f));
            obj->shader->setFloat("refractionRatio", 0.0f);
        }
        else if (type == 2) {
            obj->setShader(new Shader(vertexShaderPath, fragmentShaderPath));
            obj->shader->use();
            obj->shader->setFloat("refractionRatio", fuzzOrIOR);
            obj->shader->setFloat("averageSlope", 0.0f);
        }
        else if (type == 3) {
            obj->setShader(new Shader(vertexShaderPath, fragmentShaderPath));
        }
        else{
            // unsupported material type, use lambertian with red color
            std::cout<<"unsupported material type"<<std::endl;
            obj->setShader(new Shader(vertexShaderPath, fragmentShaderPath));
        }
        obj->shader->use();
        obj->shader->setInt("MaterialType", type);
    }
    
    // default rendering components
    std::shared_ptr<GMVPObject> obj = nullptr; // the default object, it is used for forward rendering in the default renderer
    
    int materialType; // 0: lambertian, 1: metal, 2: dielectric, 3: emissive

    // though the model matrix is stored in the object, it is also stored here for easy access
    glm::mat4 modelMatrix = glm::mat4(1.0f); // the default model matrix is identity matrix
    // ------------------------------
    // GPU ray tracing components
    glm::vec3 AA, BB; // the AABB of the object
    Texture * texture = nullptr;
    bool hasTexture() { return texture != nullptr; }
    Texture * getTexture() { return texture; }
    // ------------------------------

};

class Skybox : public RenderableSceneObject {
public:
    std::shared_ptr<GSkybox> skybox;
    Skybox(std::shared_ptr<GSkybox> _skybox) {
        this->skybox = _skybox;
        this->renderComponent = std::make_unique<ObjectRenderComponent>(_skybox, SKYBOX);
        SkyboxTexture * _skyboxTexture = _skybox->skyboxTexture;
        for (int i = 0; i < 6; i++) {
            Texture * texture = &_skyboxTexture->faces[i];
            if (texture->dataFormat != GL_UNSIGNED_BYTE) {
                std::cout<<"skybox texture data format is not GL_UNSIGNED_BYTE"<<std::endl;
                return;
            }
        }
    }

    
};


class ObjectRotationComponent : public IComponent {
public:
    ObjectRotationComponent(float _rotationSpeed = 3.0f) {
        this->rotationSpeed = _rotationSpeed;
    }
    void Tick(float deltaTime) override {
        // cast to CommonSceneObject
        if (owner == nullptr) {
            std::cout<<"[ObjectRotationComponent]: owner is nullptr"<<std::endl;
            return;
        }
        if (dynamic_cast<CommonSceneObject*>(owner)) {
            CommonSceneObject * commonSceneObject = dynamic_cast<CommonSceneObject*>(owner);
            commonSceneObject->setModelMatrix(glm::rotate(commonSceneObject->getModelMatrix(), glm::radians(rotationSpeed) * deltaTime, glm::vec3(0, 1, 0)));
        }

    }
private:
    float rotationSpeed = 3.0f;
};



class ObjectPeriodicTranslationComponent : public IComponent {
public:
    ObjectPeriodicTranslationComponent(glm::vec3 _translationSpeed = glm::vec3(0, 0, 0), glm::vec3 _translationRange = glm::vec3(0, 0, 0)) {
        this->translationSpeed = _translationSpeed;
        this->translationRange = _translationRange;
        this->currentTranslation = glm::vec3(0, 0, 0);
    }
    void Tick(float deltaTime) override {
        // cast to CommonSceneObject
        if (owner == nullptr) {
            std::cout<<"[ObjectRotationComponent]: owner is nullptr"<<std::endl;
            return;
        }
        if (dynamic_cast<CommonSceneObject*>(owner)) {
            CommonSceneObject * commonSceneObject = dynamic_cast<CommonSceneObject*>(owner);
            
            glm::vec3 scale = glm::vec3(commonSceneObject->getModelMatrix()[0][0], commonSceneObject->getModelMatrix()[1][1], commonSceneObject->getModelMatrix()[2][2]);
            glm::vec3 translation = translationSpeed * deltaTime;
            translation = translation / scale; // adjust the translation by the scale of the object
            if (glm::length(currentTranslation + translation) > glm::length(translationRange / scale)) {
                translation = -translation;
                translationSpeed = -translationSpeed;
            }
            commonSceneObject->setModelMatrix(glm::translate(commonSceneObject->getModelMatrix(), translation));
            currentTranslation = currentTranslation + translation;
        }
        else if (dynamic_cast<PointLight*>(owner)) {
            PointLight * pointLight = dynamic_cast<PointLight*>(owner);
            glm::vec3 translation = translationSpeed * deltaTime;
            if (glm::length(currentTranslation + translation) > glm::length(translationRange)) {
                translation = -translation;
                translationSpeed = -translationSpeed;
            }
            pointLight->position = pointLight->position + translation;
            currentTranslation = currentTranslation + translation;
        
        }

    }
private:
    glm::vec3 translationSpeed = glm::vec3(0, 0, 0);
    glm::vec3 translationRange = glm::vec3(0, 0, 0);
    glm::vec3 currentTranslation = glm::vec3(0, 0, 0);
};






#endif