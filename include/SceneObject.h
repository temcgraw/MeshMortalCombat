#ifndef SCENE_OBJECT_H
#define SCENE_OBJECT_H
#include <vector>
#include <Component.h>
#include <memory>


class SceneObject {
public:
    virtual ~SceneObject() {}

    void addComponent(std::unique_ptr<IComponent> component) {
        // each component should only belong to one object
        // so we should check if the component is already attached to another object
        if (component->hasOwner()) {
            std::cout << "ERROR: [addComponent] Component already has an owner\n";
            return;
        }
        component->setOwner(this);
        components.push_back(std::move(component));
    }

    virtual void Tick(float deltaTime) {
        if (!active) return;
        for (auto& component : components) {
            component->Tick(deltaTime);
        }
    }

    virtual void Active() {
        this->active = true;
        for (auto& component : components) {
            component->active = true;
        }
    }

    virtual void Deactive() {
        this->active = false;
        for (auto& component : components) {
            component->active = false;
        }
    }

    bool isActive() {
        return active;
    }

protected:
    std::vector<std::unique_ptr<IComponent>> components;
    bool active = true;
};

enum LightType {
        POINT_LIGHT,
        DIRECTIONAL_LIGHT,
        SPOT_LIGHT
};

class Light : public SceneObject {
    public:
    Light(glm::vec3 _color) : color(_color) {}
    Light() {}
    virtual ~Light() {}
    glm::vec3 color;
    enum LightType lightType;
};

class PointLight : public Light {
public:
    PointLight(glm::vec3 _position, glm::vec3 _color) : position(_position) {
        lightType = POINT_LIGHT;
        color = _color;
    }
    void set(glm::vec3 _position, glm::vec3 _color) {// update the light
        position = _position;
        color = _color;
        lightType = POINT_LIGHT;
    }
    glm::vec3 position;
};



class RenderableSceneObject : public virtual SceneObject {
public:

    virtual void attachToSceneRenderList(std::vector<std::shared_ptr<RenderComponent>>& renderQueue) {
        if (this->renderComponent == nullptr) {
            std::cout<< "ERROR: [attach] RenderComponent is nullptr\n"<<std::endl;
            return;
        }
        renderQueue.push_back(this->renderComponent);
    }
    virtual void detachFromSceneRenderList(std::vector<std::shared_ptr<RenderComponent>>& renderQueue) {
        if (this->renderComponent == nullptr) {
            std::cout<< "ERROR: [detach] RenderComponent is nullptr\n"<<std::endl;
            return;
        }
        auto it = std::find(renderQueue.begin(), renderQueue.end(), this->renderComponent);
        if (it != renderQueue.end()) {
            renderQueue.erase(it);
        }
    }

    virtual void Active() override {
        SceneObject::Active();
        this->renderComponent->active = true;
    }

    virtual void Deactive() override {
        SceneObject::Deactive();
        this->renderComponent->active = false;
    }

protected:
    std::shared_ptr<RenderComponent> renderComponent;
};


// since compute shader objects are usually highly specialized, they should be implemented separately
// so here we only provide a base class for compute and offers interfaces for the derived classes
class ComputeSceneObject : public virtual SceneObject {
public:
    virtual void attachToSceneComputeList(std::vector<std::shared_ptr<ComputeComponent>>& computeQueue) = 0;
    virtual void detachFromSceneComputeList(std::vector<std::shared_ptr<ComputeComponent>>& computeQueue) = 0;
};














#endif