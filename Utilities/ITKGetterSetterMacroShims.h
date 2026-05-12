/*=========================================================================
 *
 *  ITKGetterSetterMacroShims.h
 *
 *  Backward-compatibility shim providing the itkVirtual* / itkFinal* /
 *  itkNonVirtual* Set/Get macro family from
 *  InsightSoftwareConsortium/ITK PR #6253.
 *
 *  Once ANTs's minimum ITK version contains PR #6253, the
 *  upstream itkMacro.h itself defines itkSetMacroImpl (and the rest of
 *  the *Impl helpers plus their public variants). The sentinel guard
 *  below collapses this file to an empty translation unit in that
 *  configuration so it can be retired without source churn.
 *
 *  Usage: include AFTER itkMacro.h.
 *
 *=========================================================================*/
#ifndef ITKGetterSetterMacroShims_h
#define ITKGetterSetterMacroShims_h

#include "itkMacro.h"

// Sentinel: itkSetMacroImpl is the first internal helper introduced by
// ITK PR #6253. When ITK provides it, the entire shim is a no-op.
#if !defined(itkSetMacroImpl)

// ---------------------------------------------------------------------------
// itkSetInputMacro family
// ---------------------------------------------------------------------------
#define itkSetInputMacroImpl(virtualKeyword, finalKeyword, name, type)                   \
  virtualKeyword void Set##name(const type * _arg) finalKeyword                          \
  {                                                                                      \
    itkDebugMacro("setting input " #name " to " << _arg);                                \
    if (_arg != itkDynamicCastInDebugMode<type *>(this->ProcessObject::GetInput(#name))) \
    {                                                                                    \
      this->ProcessObject::SetInput(#name, const_cast<type *>(_arg));                    \
      this->Modified();                                                                  \
    }                                                                                    \
  }                                                                                      \
  ITK_MACROEND_NOOP_STATEMENT

#define itkVirtualSetInputMacro(name, type)    itkSetInputMacroImpl(virtual, , name, type)
#define itkFinalSetInputMacro(name, type)      itkSetInputMacroImpl(, final, name, type)
#define itkNonVirtualSetInputMacro(name, type) itkSetInputMacroImpl(, , name, type)

// ---------------------------------------------------------------------------
// itkGetInputMacro family
// ---------------------------------------------------------------------------
#define itkGetInputMacroImpl(virtualKeyword, finalKeyword, name, type)                         \
  virtualKeyword const type * Get##name() const finalKeyword                                   \
  {                                                                                            \
    itkDebugMacro("returning input " << #name " of " << this->ProcessObject::GetInput(#name)); \
    return itkDynamicCastInDebugMode<const type *>(this->ProcessObject::GetInput(#name));      \
  }                                                                                            \
  ITK_MACROEND_NOOP_STATEMENT

#define itkVirtualGetInputMacro(name, type)    itkGetInputMacroImpl(virtual, , name, type)
#define itkFinalGetInputMacro(name, type)      itkGetInputMacroImpl(, final, name, type)
#define itkNonVirtualGetInputMacro(name, type) itkGetInputMacroImpl(, , name, type)

// ---------------------------------------------------------------------------
// itkSetDecoratedInputMacro family
// ---------------------------------------------------------------------------
// clang-format off
#define itkSetDecoratedInputMacroImpl(virtualKeyword, finalKeyword, name, type)                                     \
  virtualKeyword void Set##name##Input(const SimpleDataObjectDecorator<type> * _arg) finalKeyword                   \
  {                                                                                                                 \
    itkDebugMacro("setting input " #name " to " << _arg);                                                           \
    if (_arg != itkDynamicCastInDebugMode<SimpleDataObjectDecorator<type> *>(this->ProcessObject::GetInput(#name))) \
    {                                                                                                               \
      this->ProcessObject::SetInput(#name, const_cast<SimpleDataObjectDecorator<type> *>(_arg));                    \
      this->Modified();                                                                                             \
    }                                                                                                               \
  }                                                                                                                 \
  virtualKeyword void Set##name(const SimpleDataObjectDecorator<type> * _arg) finalKeyword                          \
  {                                                                                                                 \
    this->Set##name##Input(_arg);                                                                                   \
  }                                                                                                                 \
  virtualKeyword void Set##name(const type & _arg) finalKeyword                                                     \
  {                                                                                                                 \
    using DecoratorType = SimpleDataObjectDecorator<type>;                                                          \
    itkDebugMacro("setting input " #name " to " << _arg);                                                           \
    const DecoratorType * oldInput =                                                                                \
      itkDynamicCastInDebugMode<const DecoratorType *>(this->ProcessObject::GetInput(#name));                       \
    ITK_GCC_PRAGMA_PUSH                                                                                             \
    ITK_GCC_SUPPRESS_Wfloat_equal                                                                                   \
    if (oldInput && oldInput->Get() == _arg)                                                                        \
    {                                                                                                               \
      return;                                                                                                       \
    }                                                                                                               \
    ITK_GCC_PRAGMA_POP                                                                                              \
    auto newInput = DecoratorType::New();                                                                           \
    newInput->Set(_arg);                                                                                            \
    this->Set##name##Input(newInput);                                                                               \
  }                                                                                                                 \
  ITK_MACROEND_NOOP_STATEMENT
// clang-format on

#define itkVirtualSetDecoratedInputMacro(name, type)    itkSetDecoratedInputMacroImpl(virtual, , name, type)
#define itkFinalSetDecoratedInputMacro(name, type)      itkSetDecoratedInputMacroImpl(, final, name, type)
#define itkNonVirtualSetDecoratedInputMacro(name, type) itkSetDecoratedInputMacroImpl(, , name, type)

// ---------------------------------------------------------------------------
// itkGetDecoratedInputMacro family
// ---------------------------------------------------------------------------
#define itkGetDecoratedInputMacroImpl(virtualKeyword, finalKeyword, name, type)                                      \
  virtualKeyword const SimpleDataObjectDecorator<type> * Get##name##Input() const finalKeyword                       \
  {                                                                                                                  \
    itkDebugMacro("returning input " << #name " of " << this->ProcessObject::GetInput(#name));                       \
    return itkDynamicCastInDebugMode<const SimpleDataObjectDecorator<type> *>(this->ProcessObject::GetInput(#name)); \
  }                                                                                                                  \
  virtualKeyword const type & Get##name() const finalKeyword                                                         \
  {                                                                                                                  \
    itkDebugMacro("Getting input " #name);                                                                           \
    using DecoratorType = SimpleDataObjectDecorator<type>;                                                           \
    const DecoratorType * input =                                                                                    \
      itkDynamicCastInDebugMode<const DecoratorType *>(this->ProcessObject::GetInput(#name));                        \
    if (input == nullptr)                                                                                            \
    {                                                                                                                \
      itkExceptionMacro("input" #name " is not set");                                                                \
    }                                                                                                                \
    return input->Get();                                                                                             \
  }                                                                                                                  \
  ITK_MACROEND_NOOP_STATEMENT

#define itkVirtualGetDecoratedInputMacro(name, type)    itkGetDecoratedInputMacroImpl(virtual, , name, type)
#define itkFinalGetDecoratedInputMacro(name, type)      itkGetDecoratedInputMacroImpl(, final, name, type)
#define itkNonVirtualGetDecoratedInputMacro(name, type) itkGetDecoratedInputMacroImpl(, , name, type)

// Combined Set/Get for a decorated input.
#define itkVirtualSetGetDecoratedInputMacro(name, type) \
  itkVirtualSetDecoratedInputMacro(name, type);         \
  itkVirtualGetDecoratedInputMacro(name, type)
#define itkFinalSetGetDecoratedInputMacro(name, type) \
  itkFinalSetDecoratedInputMacro(name, type);         \
  itkFinalGetDecoratedInputMacro(name, type)
#define itkNonVirtualSetGetDecoratedInputMacro(name, type) \
  itkNonVirtualSetDecoratedInputMacro(name, type);         \
  itkNonVirtualGetDecoratedInputMacro(name, type)

// ---------------------------------------------------------------------------
// itkSetDecoratedObjectInputMacro family
// ---------------------------------------------------------------------------
#define itkSetDecoratedObjectInputMacroImpl(virtualKeyword, finalKeyword, name, type)                         \
  virtualKeyword void Set##name##Input(const DataObjectDecorator<type> * _arg) finalKeyword                   \
  {                                                                                                           \
    itkDebugMacro("setting input " #name " to " << _arg);                                                     \
    if (_arg != itkDynamicCastInDebugMode<DataObjectDecorator<type> *>(this->ProcessObject::GetInput(#name))) \
    {                                                                                                         \
      this->ProcessObject::SetInput(#name, const_cast<DataObjectDecorator<type> *>(_arg));                    \
      this->Modified();                                                                                       \
    }                                                                                                         \
  }                                                                                                           \
  virtualKeyword void Set##name(const type * _arg) finalKeyword                                               \
  {                                                                                                           \
    using DecoratorType = DataObjectDecorator<type>;                                                          \
    itkDebugMacro("setting input " #name " to " << _arg);                                                     \
    const DecoratorType * oldInput =                                                                          \
      itkDynamicCastInDebugMode<const DecoratorType *>(this->ProcessObject::GetInput(#name));                 \
    if (oldInput && oldInput->Get() == _arg)                                                                  \
    {                                                                                                         \
      return;                                                                                                 \
    }                                                                                                         \
    auto newInput = DecoratorType::New();                                                                     \
    newInput->Set(_arg);                                                                                      \
    this->Set##name##Input(newInput);                                                                         \
  }                                                                                                           \
  ITK_MACROEND_NOOP_STATEMENT

#define itkVirtualSetDecoratedObjectInputMacro(name, type)    itkSetDecoratedObjectInputMacroImpl(virtual, , name, type)
#define itkFinalSetDecoratedObjectInputMacro(name, type)      itkSetDecoratedObjectInputMacroImpl(, final, name, type)
#define itkNonVirtualSetDecoratedObjectInputMacro(name, type) itkSetDecoratedObjectInputMacroImpl(, , name, type)

// ---------------------------------------------------------------------------
// itkGetDecoratedObjectInputMacro family
// ---------------------------------------------------------------------------
#define itkGetDecoratedObjectInputMacroImpl(virtualKeyword, finalKeyword, name, type)                          \
  virtualKeyword const DataObjectDecorator<type> * Get##name##Input() const finalKeyword                       \
  {                                                                                                            \
    itkDebugMacro("returning input " << #name " of " << this->ProcessObject::GetInput(#name));                 \
    return itkDynamicCastInDebugMode<const DataObjectDecorator<type> *>(this->ProcessObject::GetInput(#name)); \
  }                                                                                                            \
  virtualKeyword const type * Get##name() const finalKeyword                                                   \
  {                                                                                                            \
    itkDebugMacro("Getting input " #name);                                                                     \
    using DecoratorType = DataObjectDecorator<type>;                                                           \
    const DecoratorType * input =                                                                              \
      itkDynamicCastInDebugMode<const DecoratorType *>(this->ProcessObject::GetInput(#name));                  \
    if (input == nullptr)                                                                                      \
    {                                                                                                          \
      return nullptr;                                                                                          \
    }                                                                                                          \
    return input->Get();                                                                                       \
  }                                                                                                            \
  ITK_MACROEND_NOOP_STATEMENT

#define itkVirtualGetDecoratedObjectInputMacro(name, type)    itkGetDecoratedObjectInputMacroImpl(virtual, , name, type)
#define itkFinalGetDecoratedObjectInputMacro(name, type)      itkGetDecoratedObjectInputMacroImpl(, final, name, type)
#define itkNonVirtualGetDecoratedObjectInputMacro(name, type) itkGetDecoratedObjectInputMacroImpl(, , name, type)

#define itkVirtualSetGetDecoratedObjectInputMacro(name, type) \
  itkVirtualSetDecoratedObjectInputMacro(name, type);         \
  itkVirtualGetDecoratedObjectInputMacro(name, type)
#define itkFinalSetGetDecoratedObjectInputMacro(name, type) \
  itkFinalSetDecoratedObjectInputMacro(name, type);         \
  itkFinalGetDecoratedObjectInputMacro(name, type)
#define itkNonVirtualSetGetDecoratedObjectInputMacro(name, type) \
  itkNonVirtualSetDecoratedObjectInputMacro(name, type);         \
  itkNonVirtualGetDecoratedObjectInputMacro(name, type)

// ---------------------------------------------------------------------------
// itkSetMacro family
// ---------------------------------------------------------------------------
// clang-format off
#define itkSetMacroImpl(virtualKeyword, finalKeyword, name, type) \
  virtualKeyword void Set##name(type _arg) finalKeyword           \
  {                                                               \
    itkDebugMacro("setting " #name " to " << _arg);               \
    ITK_GCC_PRAGMA_PUSH                                           \
    ITK_GCC_SUPPRESS_Wfloat_equal                                 \
    if (this->m_##name != _arg)                                   \
    {                                                             \
      this->m_##name = std::move(_arg);                           \
      this->Modified();                                           \
    }                                                             \
    ITK_GCC_PRAGMA_POP                                            \
  }                                                               \
  ITK_MACROEND_NOOP_STATEMENT
// clang-format on

#define itkVirtualSetMacro(name, type)    itkSetMacroImpl(virtual, , name, type)
#define itkFinalSetMacro(name, type)      itkSetMacroImpl(, final, name, type)
#define itkNonVirtualSetMacro(name, type) itkSetMacroImpl(, , name, type)

// ---------------------------------------------------------------------------
// itkGetMacro / itkGetConstMacro / itkGetConstReferenceMacro families
// ---------------------------------------------------------------------------
#define itkGetMacroImpl(virtualKeyword, finalKeyword, name, type)         \
  virtualKeyword type Get##name() finalKeyword { return this->m_##name; } \
  ITK_MACROEND_NOOP_STATEMENT

#define itkGetConstMacroImpl(virtualKeyword, finalKeyword, name, type)          \
  virtualKeyword type Get##name() const finalKeyword { return this->m_##name; } \
  ITK_MACROEND_NOOP_STATEMENT

#define itkGetConstReferenceMacroImpl(virtualKeyword, finalKeyword, name, type)         \
  virtualKeyword const type & Get##name() const finalKeyword { return this->m_##name; } \
  ITK_MACROEND_NOOP_STATEMENT

#define itkVirtualGetMacro(name, type)    itkGetMacroImpl(virtual, , name, type)
#define itkFinalGetMacro(name, type)      itkGetMacroImpl(, final, name, type)
#define itkNonVirtualGetMacro(name, type) itkGetMacroImpl(, , name, type)

#define itkVirtualGetConstMacro(name, type)    itkGetConstMacroImpl(virtual, , name, type)
#define itkFinalGetConstMacro(name, type)      itkGetConstMacroImpl(, final, name, type)
#define itkNonVirtualGetConstMacro(name, type) itkGetConstMacroImpl(, , name, type)

#define itkVirtualGetConstReferenceMacro(name, type)    itkGetConstReferenceMacroImpl(virtual, , name, type)
#define itkFinalGetConstReferenceMacro(name, type)      itkGetConstReferenceMacroImpl(, final, name, type)
#define itkNonVirtualGetConstReferenceMacro(name, type) itkGetConstReferenceMacroImpl(, , name, type)

// ---------------------------------------------------------------------------
// itkSetEnumMacro / itkGetEnumMacro families
// ---------------------------------------------------------------------------
#define itkSetEnumMacroImpl(virtualKeyword, finalKeyword, name, type)  \
  virtualKeyword void Set##name(const type _arg) finalKeyword          \
  {                                                                    \
    itkDebugMacro("setting " #name " to " << static_cast<long>(_arg)); \
    if (this->m_##name != _arg)                                        \
    {                                                                  \
      this->m_##name = _arg;                                           \
      this->Modified();                                                \
    }                                                                  \
  }                                                                    \
  ITK_MACROEND_NOOP_STATEMENT

#define itkVirtualSetEnumMacro(name, type)    itkSetEnumMacroImpl(virtual, , name, type)
#define itkFinalSetEnumMacro(name, type)      itkSetEnumMacroImpl(, final, name, type)
#define itkNonVirtualSetEnumMacro(name, type) itkSetEnumMacroImpl(, , name, type)

#define itkGetEnumMacroImpl(virtualKeyword, finalKeyword, name, type)           \
  virtualKeyword type Get##name() const finalKeyword { return this->m_##name; } \
  ITK_MACROEND_NOOP_STATEMENT

#define itkVirtualGetEnumMacro(name, type)    itkGetEnumMacroImpl(virtual, , name, type)
#define itkFinalGetEnumMacro(name, type)      itkGetEnumMacroImpl(, final, name, type)
#define itkNonVirtualGetEnumMacro(name, type) itkGetEnumMacroImpl(, , name, type)

// ---------------------------------------------------------------------------
// itkSetStringMacro / itkGetStringMacro families
// ---------------------------------------------------------------------------
#define itkSetStringMacroImpl(virtualKeyword, finalKeyword, name)                                         \
  virtualKeyword void Set##name(const char * _arg) finalKeyword                                           \
  {                                                                                                       \
    if (_arg && (_arg == this->m_##name))                                                                 \
    {                                                                                                     \
      return;                                                                                             \
    }                                                                                                     \
    if (_arg)                                                                                             \
    {                                                                                                     \
      this->m_##name = _arg;                                                                              \
    }                                                                                                     \
    else                                                                                                  \
    {                                                                                                     \
      this->m_##name = "";                                                                                \
    }                                                                                                     \
    this->Modified();                                                                                     \
  }                                                                                                       \
  virtualKeyword void Set##name(const std::string & _arg) finalKeyword { this->Set##name(_arg.c_str()); } \
  ITK_MACROEND_NOOP_STATEMENT

#define itkVirtualSetStringMacro(name)    itkSetStringMacroImpl(virtual, , name)
#define itkFinalSetStringMacro(name)      itkSetStringMacroImpl(, final, name)
#define itkNonVirtualSetStringMacro(name) itkSetStringMacroImpl(, , name)

#define itkGetStringMacroImpl(virtualKeyword, finalKeyword, name)                               \
  virtualKeyword const char * Get##name() const finalKeyword { return this->m_##name.c_str(); } \
  ITK_MACROEND_NOOP_STATEMENT

#define itkVirtualGetStringMacro(name)    itkGetStringMacroImpl(virtual, , name)
#define itkFinalGetStringMacro(name)      itkGetStringMacroImpl(, final, name)
#define itkNonVirtualGetStringMacro(name) itkGetStringMacroImpl(, , name)

// ---------------------------------------------------------------------------
// itkSetClampMacro family
// ---------------------------------------------------------------------------
// clang-format off
#define itkSetClampMacroImpl(virtualKeyword, finalKeyword, name, type, min, max) \
  virtualKeyword void Set##name(type _arg) finalKeyword                          \
  {                                                                              \
    const type temp_extrema = (_arg <= min ? min : (_arg >= max ? max : _arg));  \
    itkDebugMacro("setting " << #name " to " << _arg);                           \
    ITK_GCC_PRAGMA_PUSH                                                          \
    ITK_GCC_SUPPRESS_Wfloat_equal                                                \
    if (this->m_##name != temp_extrema)                                          \
    {                                                                            \
      this->m_##name = temp_extrema;                                             \
      this->Modified();                                                          \
    }                                                                            \
    ITK_GCC_PRAGMA_POP                                                           \
  }                                                                              \
  ITK_MACROEND_NOOP_STATEMENT
// clang-format on

#define itkVirtualSetClampMacro(name, type, min, max)    itkSetClampMacroImpl(virtual, , name, type, min, max)
#define itkFinalSetClampMacro(name, type, min, max)      itkSetClampMacroImpl(, final, name, type, min, max)
#define itkNonVirtualSetClampMacro(name, type, min, max) itkSetClampMacroImpl(, , name, type, min, max)

// ---------------------------------------------------------------------------
// itkSetObjectMacro / itkGetObjectMacro / itkGetConstObjectMacro /
// itkGetModifiableObjectMacro / itkGetConstReferenceObjectMacro /
// itkSetConstObjectMacro families
// ---------------------------------------------------------------------------
#define itkSetObjectMacroImpl(virtualKeyword, finalKeyword, name, type) \
  virtualKeyword void Set##name(type * _arg) finalKeyword               \
  {                                                                     \
    itkDebugMacro("setting " << #name " to " << _arg);                  \
    if (this->m_##name != _arg)                                         \
    {                                                                   \
      this->m_##name = _arg;                                            \
      this->Modified();                                                 \
    }                                                                   \
  }                                                                     \
  ITK_MACROEND_NOOP_STATEMENT

#define itkVirtualSetObjectMacro(name, type)    itkSetObjectMacroImpl(virtual, , name, type)
#define itkFinalSetObjectMacro(name, type)      itkSetObjectMacroImpl(, final, name, type)
#define itkNonVirtualSetObjectMacro(name, type) itkSetObjectMacroImpl(, , name, type)

#define itkGetConstObjectMacroImpl(virtualKeyword, finalKeyword, name, type)                         \
  virtualKeyword const type * Get##name() const finalKeyword { return this->m_##name.GetPointer(); } \
  ITK_MACROEND_NOOP_STATEMENT

#define itkVirtualGetConstObjectMacro(name, type)    itkGetConstObjectMacroImpl(virtual, , name, type)
#define itkFinalGetConstObjectMacro(name, type)      itkGetConstObjectMacroImpl(, final, name, type)
#define itkNonVirtualGetConstObjectMacro(name, type) itkGetConstObjectMacroImpl(, , name, type)

#if defined(ITK_FUTURE_LEGACY_REMOVE)
#  define itkGetObjectMacroImpl(virtualKeyword, finalKeyword, name, type)                         \
    virtualKeyword type * Get##name() finalKeyword                                                \
    {                                                                                             \
      purposeful_error("itkGetObjectMacro should be replaced with itkGetModifiableObjectMacro."); \
    }                                                                                             \
    ITK_MACROEND_NOOP_STATEMENT

#  define itkGetModifiableObjectMacroImpl(virtualKeyword, finalKeyword, name, type)                  \
    virtualKeyword type * GetModifiable##name() finalKeyword { return this->m_##name.GetPointer(); } \
    itkGetConstObjectMacroImpl(virtualKeyword, finalKeyword, name, type)
#else
#  define itkGetObjectMacroImpl(virtualKeyword, finalKeyword, name, type)                  \
    virtualKeyword type * Get##name() finalKeyword { return this->m_##name.GetPointer(); } \
    ITK_MACROEND_NOOP_STATEMENT

#  define itkGetModifiableObjectMacroImpl(virtualKeyword, finalKeyword, name, type)                  \
    virtualKeyword type * GetModifiable##name() finalKeyword { return this->m_##name.GetPointer(); } \
    itkGetConstObjectMacroImpl(virtualKeyword, finalKeyword, name, type);                            \
    itkGetObjectMacroImpl(virtualKeyword, finalKeyword, name, type)
#endif

#define itkVirtualGetObjectMacro(name, type)    itkGetObjectMacroImpl(virtual, , name, type)
#define itkFinalGetObjectMacro(name, type)      itkGetObjectMacroImpl(, final, name, type)
#define itkNonVirtualGetObjectMacro(name, type) itkGetObjectMacroImpl(, , name, type)

#define itkVirtualGetModifiableObjectMacro(name, type)    itkGetModifiableObjectMacroImpl(virtual, , name, type)
#define itkFinalGetModifiableObjectMacro(name, type)      itkGetModifiableObjectMacroImpl(, final, name, type)
#define itkNonVirtualGetModifiableObjectMacro(name, type) itkGetModifiableObjectMacroImpl(, , name, type)

#define itkGetConstReferenceObjectMacroImpl(virtualKeyword, finalKeyword, name, type)                     \
  virtualKeyword const typename type::Pointer & Get##name() const finalKeyword { return this->m_##name; } \
  ITK_MACROEND_NOOP_STATEMENT

#define itkVirtualGetConstReferenceObjectMacro(name, type)    itkGetConstReferenceObjectMacroImpl(virtual, , name, type)
#define itkFinalGetConstReferenceObjectMacro(name, type)      itkGetConstReferenceObjectMacroImpl(, final, name, type)
#define itkNonVirtualGetConstReferenceObjectMacro(name, type) itkGetConstReferenceObjectMacroImpl(, , name, type)

#define itkSetConstObjectMacroImpl(virtualKeyword, finalKeyword, name, type) \
  virtualKeyword void Set##name(const type * _arg) finalKeyword              \
  {                                                                          \
    itkDebugMacro("setting " << #name " to " << _arg);                       \
    if (this->m_##name != _arg)                                              \
    {                                                                        \
      this->m_##name = _arg;                                                 \
      this->Modified();                                                      \
    }                                                                        \
  }                                                                          \
  ITK_MACROEND_NOOP_STATEMENT

#define itkVirtualSetConstObjectMacro(name, type)    itkSetConstObjectMacroImpl(virtual, , name, type)
#define itkFinalSetConstObjectMacro(name, type)      itkSetConstObjectMacroImpl(, final, name, type)
#define itkNonVirtualSetConstObjectMacro(name, type) itkSetConstObjectMacroImpl(, , name, type)

// ---------------------------------------------------------------------------
// itkBooleanMacro family
// ---------------------------------------------------------------------------
#define itkBooleanMacroImpl(virtualKeyword, finalKeyword, name)            \
  virtualKeyword void name##On() finalKeyword { this->Set##name(true); }   \
  virtualKeyword void name##Off() finalKeyword { this->Set##name(false); } \
  ITK_MACROEND_NOOP_STATEMENT

#define itkVirtualBooleanMacro(name)    itkBooleanMacroImpl(virtual, , name)
#define itkFinalBooleanMacro(name)      itkBooleanMacroImpl(, final, name)
#define itkNonVirtualBooleanMacro(name) itkBooleanMacroImpl(, , name)

// ---------------------------------------------------------------------------
// itkSetVectorMacro / itkGetVectorMacro families
// ---------------------------------------------------------------------------
#define itkSetVectorMacroImpl(virtualKeyword, finalKeyword, name, type, count) \
  virtualKeyword void Set##name(type data[]) finalKeyword                      \
  {                                                                            \
    if (ContainerCopyWithCheck(this->m_##name, data, count))                   \
    {                                                                          \
      this->Modified();                                                        \
    }                                                                          \
  }                                                                            \
  ITK_MACROEND_NOOP_STATEMENT

#define itkVirtualSetVectorMacro(name, type, count)    itkSetVectorMacroImpl(virtual, , name, type, count)
#define itkFinalSetVectorMacro(name, type, count)      itkSetVectorMacroImpl(, final, name, type, count)
#define itkNonVirtualSetVectorMacro(name, type, count) itkSetVectorMacroImpl(, , name, type, count)

#define itkGetVectorMacroImpl(virtualKeyword, finalKeyword, name, type, count)    \
  virtualKeyword type * Get##name() const finalKeyword { return this->m_##name; } \
  ITK_MACROEND_NOOP_STATEMENT

#define itkVirtualGetVectorMacro(name, type, count)    itkGetVectorMacroImpl(virtual, , name, type, count)
#define itkFinalGetVectorMacro(name, type, count)      itkGetVectorMacroImpl(, final, name, type, count)
#define itkNonVirtualGetVectorMacro(name, type, count) itkGetVectorMacroImpl(, , name, type, count)

// ---------------------------------------------------------------------------
// itkSetDecoratedOutputMacro / itkGetDecoratedOutputMacro families
// ---------------------------------------------------------------------------
#define itkSetDecoratedOutputMacroImpl(virtualKeyword, finalKeyword, name, type)                                     \
  virtualKeyword void Set##name##Output(const SimpleDataObjectDecorator<type> * _arg) finalKeyword                   \
  {                                                                                                                  \
    itkDebugMacro("setting output " #name " to " << _arg);                                                           \
    if (_arg != itkDynamicCastInDebugMode<SimpleDataObjectDecorator<type> *>(this->ProcessObject::GetOutput(#name))) \
    {                                                                                                                \
      this->ProcessObject::SetOutput(#name, const_cast<SimpleDataObjectDecorator<type> *>(_arg));                    \
      this->Modified();                                                                                              \
    }                                                                                                                \
  }                                                                                                                  \
  virtualKeyword void Set##name(const type & _arg) finalKeyword                                                      \
  {                                                                                                                  \
    using DecoratorType = SimpleDataObjectDecorator<type>;                                                           \
    itkDebugMacro("setting output " #name " to " << _arg);                                                           \
    DecoratorType * output =                                                                                         \
      itkDynamicCastInDebugMode<DecoratorType *>(this->ProcessObject::GetOutput(#name));                             \
    if (output)                                                                                                      \
    {                                                                                                                \
      if (output->Get() == _arg)                                                                                     \
      {                                                                                                              \
        return;                                                                                                      \
      }                                                                                                              \
      else                                                                                                           \
      {                                                                                                              \
        output->Set(_arg);                                                                                           \
      }                                                                                                              \
    }                                                                                                                \
    else                                                                                                             \
    {                                                                                                                \
      auto newOutput = DecoratorType::New();                                                                         \
      newOutput->Set(_arg);                                                                                          \
      this->Set##name##Output(newOutput);                                                                            \
    }                                                                                                                \
  }                                                                                                                  \
  ITK_MACROEND_NOOP_STATEMENT

#define itkVirtualSetDecoratedOutputMacro(name, type)    itkSetDecoratedOutputMacroImpl(virtual, , name, type)
#define itkFinalSetDecoratedOutputMacro(name, type)      itkSetDecoratedOutputMacroImpl(, final, name, type)
#define itkNonVirtualSetDecoratedOutputMacro(name, type) itkSetDecoratedOutputMacroImpl(, , name, type)

#define itkGetDecoratedOutputMacroImpl(virtualKeyword, finalKeyword, name, type)                                      \
  virtualKeyword const SimpleDataObjectDecorator<type> * Get##name##Output() const finalKeyword                       \
  {                                                                                                                   \
    itkDebugMacro("returning output " << #name " of " << this->ProcessObject::GetOutput(#name));                      \
    return itkDynamicCastInDebugMode<const SimpleDataObjectDecorator<type> *>(this->ProcessObject::GetOutput(#name)); \
  }                                                                                                                   \
  virtualKeyword const type & Get##name() const finalKeyword                                                          \
  {                                                                                                                   \
    itkDebugMacro("Getting output " #name);                                                                           \
    using DecoratorType = SimpleDataObjectDecorator<type>;                                                            \
    const DecoratorType * output =                                                                                    \
      itkDynamicCastInDebugMode<const DecoratorType *>(this->ProcessObject::GetOutput(#name));                        \
    if (output == nullptr)                                                                                            \
    {                                                                                                                 \
      itkExceptionMacro("output" #name " is not set");                                                                \
    }                                                                                                                 \
    return output->Get();                                                                                             \
  }                                                                                                                   \
  ITK_MACROEND_NOOP_STATEMENT

#define itkVirtualGetDecoratedOutputMacro(name, type)    itkGetDecoratedOutputMacroImpl(virtual, , name, type)
#define itkFinalGetDecoratedOutputMacro(name, type)      itkGetDecoratedOutputMacroImpl(, final, name, type)
#define itkNonVirtualGetDecoratedOutputMacro(name, type) itkGetDecoratedOutputMacroImpl(, , name, type)

#endif // !defined(itkSetMacroImpl)

#endif // ITKGetterSetterMacroShims_h
