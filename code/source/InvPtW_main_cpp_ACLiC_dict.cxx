// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIUsersdIstephanstiefelmaierdIDocumentsdIPromotiondIrepos_clouddIprojectsdI2024dI2024mI03mI06_investigatePtWeightsdIdOdIsourcedIInvPtW_main_cpp_ACLiC_dict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/Users/stephanstiefelmaier/Documents/Promotion/repos_cloud/projects/2024/2024-03-06_investigatePtWeights/./source/InvPtW_main.cpp"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *InvPtW_main_Dictionary();
   static void InvPtW_main_TClassManip(TClass*);
   static void delete_InvPtW_main(void *p);
   static void deleteArray_InvPtW_main(void *p);
   static void destruct_InvPtW_main(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::InvPtW_main*)
   {
      ::InvPtW_main *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::InvPtW_main));
      static ::ROOT::TGenericClassInfo 
         instance("InvPtW_main", "source/InvPtW_main.h", 26,
                  typeid(::InvPtW_main), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &InvPtW_main_Dictionary, isa_proxy, 4,
                  sizeof(::InvPtW_main) );
      instance.SetDelete(&delete_InvPtW_main);
      instance.SetDeleteArray(&deleteArray_InvPtW_main);
      instance.SetDestructor(&destruct_InvPtW_main);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::InvPtW_main*)
   {
      return GenerateInitInstanceLocal(static_cast<::InvPtW_main*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::InvPtW_main*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *InvPtW_main_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::InvPtW_main*>(nullptr))->GetClass();
      InvPtW_main_TClassManip(theClass);
   return theClass;
   }

   static void InvPtW_main_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","/Users/stephanstiefelmaier/Documents/Promotion/repos_cloud/projects/2024/2024-03-06_investigatePtWeights/./source/InvPtW_main.h");
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrapper around operator delete
   static void delete_InvPtW_main(void *p) {
      delete (static_cast<::InvPtW_main*>(p));
   }
   static void deleteArray_InvPtW_main(void *p) {
      delete [] (static_cast<::InvPtW_main*>(p));
   }
   static void destruct_InvPtW_main(void *p) {
      typedef ::InvPtW_main current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::InvPtW_main

namespace {
  void TriggerDictionaryInitialization_InvPtW_main_cpp_ACLiC_dict_Impl() {
    static const char* headers[] = {
"./source/InvPtW_main.cpp",
nullptr
    };
    static const char* includePaths[] = {
"/Users/stephanstiefelmaier/work/alice/sw/osx_x86-64/ROOT/v6-30-01-alice4-local1/include",
"/Users/stephanstiefelmaier/work/alice/sw/osx_x86-64/TBB/v2021.5.0-local1/include",
"/Users/stephanstiefelmaier/work/alice/sw/osx_x86-64/boost/v1.83.0-alice2-local1/include",
"/Users/stephanstiefelmaier/work/alice/sw/osx_x86-64/ROOT/v6-30-01-alice4-local1/etc/",
"/Users/stephanstiefelmaier/work/alice/sw/osx_x86-64/ROOT/v6-30-01-alice4-local1/etc//cling",
"/Users/stephanstiefelmaier/work/alice/sw/osx_x86-64/ROOT/v6-30-01-alice4-local1/etc//cling/plugins/include",
"/Users/stephanstiefelmaier/work/alice/sw/osx_x86-64/ROOT/v6-30-01-alice4-local1/include/",
"/Users/stephanstiefelmaier/work/alice/sw/osx_x86-64/ROOT/v6-30-01-alice4-local1/include",
"/Users/stephanstiefelmaier/work/alice/sw/osx_x86-64/ROOT/v6-30-01-alice4-local1/include/",
"/Users/stephanstiefelmaier/Documents/Promotion/repos_cloud/projects/2024/2024-03-06_investigatePtWeights/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "InvPtW_main_cpp_ACLiC_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$./source/InvPtW_main.cpp")))  InvPtW_main;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "InvPtW_main_cpp_ACLiC_dict dictionary payload"

#ifndef __ACLIC__
  #define __ACLIC__ 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "./source/InvPtW_main.cpp"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"", payloadCode, "@",
"InvPtW_main", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("InvPtW_main_cpp_ACLiC_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_InvPtW_main_cpp_ACLiC_dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_InvPtW_main_cpp_ACLiC_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_InvPtW_main_cpp_ACLiC_dict() {
  TriggerDictionaryInitialization_InvPtW_main_cpp_ACLiC_dict_Impl();
}
