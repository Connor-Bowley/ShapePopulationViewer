#include "ShapePopulationData.h"
#include <vtkTriangle.h>

// Relax JSON standard and allow reading/writing of nan and inf
// values. Such values should not normally occur, but if they do then
// it is easier to troubleshoot problems if numerical values are the
// same in memory and files.
// kWriteNanAndInfFlag = 2,        //!< Allow writing of Infinity, -Infinity and NaN.
#define RAPIDJSON_WRITE_DEFAULT_FLAGS 2
// kParseNanAndInfFlag = 256,      //!< Allow parsing NaN, Inf, Infinity, -Inf and -Infinity as doubles.
#define RAPIDJSON_PARSE_DEFAULT_FLAGS 256

#include <rapidjson/document.h>     // rapidjson's DOM-style API
#include <rapidjson/prettywriter.h> // for stringify JSON
#include <rapidjson/filereadstream.h>
#include <rapidjson/filewritestream.h>

static bool endswith(std::string file, std::string ext)
{
    int epos = file.length() - ext.length();
    if (epos < 0)
    {
        return false;
    }
    return file.rfind(ext) == (unsigned int)epos;
}

ShapePopulationData::ShapePopulationData()
{
    m_PolyData = vtkSmartPointer<vtkPolyData>::New();
}

vtkSmartPointer<vtkPolyData> ShapePopulationData::ReadPolyData(std::string a_filePath)
{
    if (endswith(a_filePath, ".vtp"))
    {
        vtkSmartPointer<vtkXMLPolyDataReader> meshReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
        meshReader->SetFileName(a_filePath.c_str());
        meshReader->Update();
        return meshReader->GetOutput();
    }
    else if (endswith(a_filePath, ".vtk"))
    {
        vtkSmartPointer<vtkPolyDataReader> meshReader = vtkSmartPointer<vtkPolyDataReader>::New();
        meshReader->SetFileName(a_filePath.c_str());
        meshReader->Update();
        return meshReader->GetOutput();
    }
    else
    {
        return NULL;
    }
}

vtkSmartPointer<vtkPolyData> ShapePopulationData::ReadMesh(std::string a_filePath)
{
    if (endswith(a_filePath, ".vtp") || endswith(a_filePath, ".vtk"))
    {
        vtkSmartPointer<vtkPolyData> polyData = ReadPolyData(a_filePath);
        if(polyData == NULL) return 0;

        ReadMesh(polyData, a_filePath);
    }
    else if (endswith(a_filePath, ".xml"))
    {
        vtkNew<vtkXMLDataParser> parser;
        parser->SetFileName(a_filePath.c_str());
        if( parser->Parse() != 1) return 0;

        auto* root = parser->GetRootElement();
        vtkSmartPointer<vtkPolyData> upPD = ReadPolyData(root->FindNestedElementWithName("upSpoke")->GetCharacterData());
        vtkSmartPointer<vtkPolyData> downPD = ReadPolyData(root->FindNestedElementWithName("downSpoke")->GetCharacterData());
        vtkSmartPointer<vtkPolyData> crestPD = ReadPolyData(root->FindNestedElementWithName("crestSpoke")->GetCharacterData());
        if(upPD == NULL || downPD == NULL || crestPD == NULL) return 0;

        ReadSRep(upPD, downPD, crestPD, a_filePath);
    }
    else if (endswith(a_filePath, ".srep.json"))
    {
        auto srepPolyData = ReadSRepJson(a_filePath);
        if (srepPolyData)
        {
            m_PolyData = srepPolyData;
            int numAttributes = m_PolyData->GetPointData()->GetNumberOfArrays();
            for (int j = 0; j < numAttributes; j++)
            {
                int dim = m_PolyData->GetPointData()->GetArray(j)->GetNumberOfComponents();
                const char * AttributeName = m_PolyData->GetPointData()->GetArrayName(j);

                if (dim == 1)
                {
                    std::string AttributeString = AttributeName;
                    m_AttributeList.push_back(AttributeString);
                }
            }
            std::sort(m_AttributeList.begin(),m_AttributeList.end());
        }
    }
    else
    {
        return 0;
    }
    return m_PolyData;
}

void ShapePopulationData::ReadMesh(vtkPolyData* polyData, const std::string& a_filePath)
{
    size_t found = a_filePath.find_last_of("/\\");
    if (found != std::string::npos)
    {
        m_FileDir = a_filePath.substr(0,found);
        m_FileName = a_filePath.substr(found+1);
    }
    else
    {
        m_FileName = a_filePath;
    }

    vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
    normalGenerator->SetInputData(polyData);
    normalGenerator->SplittingOff();
    normalGenerator->ComputePointNormalsOn();
    normalGenerator->ComputeCellNormalsOff();
    normalGenerator->Update();

    //Update the class members
    m_PolyData = normalGenerator->GetOutput();

    int numAttributes = m_PolyData->GetPointData()->GetNumberOfArrays();
    for (int j = 0; j < numAttributes; j++)
    {
        int dim = m_PolyData->GetPointData()->GetArray(j)->GetNumberOfComponents();
        const char * AttributeName = m_PolyData->GetPointData()->GetArrayName(j);

        if (dim == 1 || dim == 3 )
        {
            std::string AttributeString = AttributeName;
            m_AttributeList.push_back(AttributeString);
        }
        if( dim == 3)
        {
            //Vectors
            vtkPVPostFilter *  getVectors = vtkPVPostFilter::New();
            std::ostringstream strs;
            strs.str("");
            strs.clear();
            strs << AttributeName << "_mag" << std::endl;
            getVectors->DoAnyNeededConversions(m_PolyData,strs.str().c_str(),vtkDataObject::FIELD_ASSOCIATION_POINTS, AttributeName, "Magnitude");
            getVectors->Delete();
        }
    }
    std::sort(m_AttributeList.begin(),m_AttributeList.end());
}

void ShapePopulationData::ReadSRep(vtkPolyData* upPD, vtkPolyData* downPD, vtkPolyData* crestPD, const std::string& a_filePath)
{
    size_t found = a_filePath.find_last_of("/\\");
    if (found != std::string::npos)
    {
        m_FileDir = a_filePath.substr(0,found);
        m_FileName = a_filePath.substr(found+1);
    }
    else
    {
        m_FileName = a_filePath;
    }

    // Create spoke geometries and merge poly data
    vtkNew<vtkAppendPolyData> append;
    append->AddInputData(ExtractSpokes(upPD, 0));
    append->AddInputData(ExtractSpokes(downPD, 4));
    append->AddInputData(ExtractSpokes(crestPD, 2));
    append->AddInputData(AttachCellType(upPD, 1));
    append->AddInputData(AttachCellType(crestPD, 3));
    append->Update();

    // Update the class members
    m_PolyData = append->GetOutput();

    int numAttributes = m_PolyData->GetPointData()->GetNumberOfArrays();
    for (int j = 0; j < numAttributes; j++)
    {
        int dim = m_PolyData->GetPointData()->GetArray(j)->GetNumberOfComponents();
        const char * AttributeName = m_PolyData->GetPointData()->GetArrayName(j);

        if (dim == 1)
        {
            std::string AttributeString = AttributeName;
            m_AttributeList.push_back(AttributeString);
        }
    }
    std::sort(m_AttributeList.begin(),m_AttributeList.end());
}

vtkSmartPointer<vtkPolyData> ShapePopulationData::ExtractSpokes(vtkPolyData* polyData, int cellType)
{
    int nPoints = polyData->GetNumberOfPoints();
    auto* pointData = polyData->GetPointData();
    auto* arr_length = pointData->GetArray("spokeLength");
    auto* arr_dirs = pointData->GetArray("spokeDirection");

    vtkNew<vtkPoints> spokePoints;
    vtkNew<vtkCellArray> spokeLines;
    vtkNew<vtkIntArray> typeArray;
    typeArray->SetName("cellType");
    typeArray->SetNumberOfComponents(1);

    for (int i  = 0; i < nPoints; ++i)
    {
        auto* pt0 = polyData->GetPoint(i);
        auto spoke_length = arr_length->GetTuple1(i);
        auto* dir = arr_dirs->GetTuple3(i);
        double pt1[] = {pt0[0] + spoke_length * dir[0],
                        pt0[1] + spoke_length * dir[1],
                        pt0[2] + spoke_length * dir[2]};
        int id0 = spokePoints->InsertNextPoint(pt0);
        int id1 = spokePoints->InsertNextPoint(pt1);

        vtkNew<vtkLine> arrow;
        arrow->GetPointIds()->SetId(0, id0);
        arrow->GetPointIds()->SetId(1, id1);
        spokeLines->InsertNextCell(arrow);
        typeArray->InsertNextValue(cellType);
        typeArray->InsertNextValue(cellType);
    }

    vtkNew<vtkPolyData> spokePD;
    spokePD->SetPoints(spokePoints);
    spokePD->SetLines(spokeLines);
    spokePD->GetPointData()->SetActiveScalars("cellType");
    spokePD->GetPointData()->SetScalars(typeArray);
    return spokePD;
}

vtkSmartPointer<vtkPolyData> ShapePopulationData::AttachCellType(vtkPolyData *polyData, int cellType)
{
    vtkNew<vtkIntArray> outputType;
    outputType->SetName("cellType");
    outputType->SetNumberOfComponents(1);
    for (int i = 0; i < polyData->GetNumberOfPoints(); ++i)
    {
        outputType->InsertNextValue(cellType);
    }

    polyData->GetPointData()->SetActiveScalars("cellType");
    polyData->GetPointData()->SetScalars(outputType);
    return polyData;
}

namespace {

class FileRAII {
public:
    FileRAII(const std::string& a_filePath)
    {
        m_Handle = fopen(a_filePath.c_str(), "r");
    }
    ~FileRAII()
    {
        fclose(m_Handle);
    }
    FILE* GetHandle()
    {
        return m_Handle;
    }
private:
    FILE* m_Handle;
};

namespace keys {
  const char * const EllipticalSRep = "EllipticalSRep";
  const char * const CrestPoints = "CrestPoints";
  const char * const Steps = "Steps";
  const char * const Skeleton = "Skeleton";
  const char * const UpSpoke = "UpSpoke";
  const char * const DownSpoke = "DownSpoke";
  const char * const CrestSpoke = "CrestSpoke";
  const char * const Direction = "Direction";
  const char * const SkeletalPoint = "SkeletalPoint";
  const char * const Value = "Value";
  const char * const CoordinateSystem = "CoordinateSystem";
}

rapidjson::Value::MemberIterator SafeFindMember(rapidjson::Value& json, const char* name)
{
  auto iter = json.FindMember(name);
  if (iter == json.MemberEnd())
  {
    throw std::invalid_argument(std::string("Error finding json member '") + name + "'");
  }
  return iter;
}

double readDouble(rapidjson::Value& json)
{
  if (!json.IsDouble())
  {
    throw std::invalid_argument("Expected a JSON double.");
  }
  return json.GetDouble();
}

int readUint(rapidjson::Value& json) {
  if (!json.IsUint()) {
    throw std::invalid_argument("Expected a JSON uint.");
  }
  return json.GetUint();
}

std::array<double, 3> read3DArray(rapidjson::Value& json)
{
    if (!json.IsArray())
    {
        throw std::invalid_argument("Attempting to read an array that is not a json array");
    }
    auto jsonArray = json.GetArray();

    if (jsonArray.Size() != 3)
    {
        throw std::invalid_argument("Attempting to read a 3D array that doesn't have 3 dimensions");
    }

    std::array<double, 3> array;
    std::transform(jsonArray.Begin(), jsonArray.End(), array.begin(), readDouble);
    return array;
}

std::array<double, 3> read3DArrayAsRAS(rapidjson::Value& json)
{
    if (!json.IsObject()) {
        throw std::invalid_argument("Attempting to read a LPS/RAS array that is not a json object");
    }
    auto value  = read3DArray(SafeFindMember(json, keys::Value)->value);
    if (SafeFindMember(json, keys::CoordinateSystem)->value == "RAS") {
        return value;
    } else {
        return std::array<double, 3>{-value[0], -value[1], value[2]};
    }
}

vtkIdType addSpoke(vtkPolyData& polyData, rapidjson::Value& json, int cellType)
{
    auto spokePoints = polyData.GetPoints();
    auto spokeLines = polyData.GetLines();
    auto typeArray = static_cast<vtkIntArray*>(polyData.GetPointData()->GetScalars("cellType"));

    auto skeletalIter = SafeFindMember(json, keys::SkeletalPoint);
    auto directionIter = SafeFindMember(json, keys::Direction);

    const auto skeletalPoint = read3DArrayAsRAS(skeletalIter->value);
    const auto directionVector = read3DArrayAsRAS(directionIter->value);
    const std::array<double, 3> boundaryPoint = {
        skeletalPoint[0] + directionVector[0],
        skeletalPoint[1] + directionVector[1],
        skeletalPoint[2] + directionVector[2]
    };
    auto idS = spokePoints->InsertNextPoint(skeletalPoint.data());
    auto idB = spokePoints->InsertNextPoint(boundaryPoint.data());

    vtkNew<vtkLine> arrow;
    arrow->GetPointIds()->SetId(0, idS);
    arrow->GetPointIds()->SetId(1, idB);
    spokeLines->InsertNextCell(arrow);

    typeArray->InsertNextValue(cellType);
    typeArray->InsertNextValue(cellType);

    return idS;
}

void addCrestCurve(vtkPolyData& polyData, const std::vector<vtkIdType>& crestIds) {
    auto spokePoints = polyData.GetPoints();
    auto spokeLines = polyData.GetLines();
    auto typeArray = static_cast<vtkIntArray*>(polyData.GetPointData()->GetScalars("cellType"));

    // add lines for the crest, but duplicate the points so we can set different scalars
    double p[3]; // buffer to store points in
    for (size_t i = 0; i < crestIds.size(); ++i)
    {
        auto id0 = crestIds[i];
        auto id1 = crestIds[(i + 1) % crestIds.size()];

        spokePoints->GetPoint(id0, p);
        auto id0Dup = spokePoints->InsertNextPoint(p);

        spokePoints->GetPoint(id1, p);
        auto id1Dup = spokePoints->InsertNextPoint(p);

        vtkNew<vtkLine> arrow;
        arrow->GetPointIds()->SetId(0, id0Dup);
        arrow->GetPointIds()->SetId(1, id1Dup);
        spokeLines->InsertNextCell(arrow);
        typeArray->InsertNextValue(3);
        typeArray->InsertNextValue(3);
    }
}

void addMedialSheet(vtkPolyData& polyData, const std::vector<vtkIdType>& upIds, unsigned int numFoldPoints, unsigned int numStepsPlusCrest) {
    auto spokePoints = polyData.GetPoints();
    auto medialSheet = polyData.GetPolys();
    auto typeArray = static_cast<vtkIntArray*>(polyData.GetPointData()->GetScalars("cellType"));
    // add the medial sheet cells, but duplicate the points so we can set different scalars
    std::vector<std::vector<vtkIdType>> newUpIds(numFoldPoints, std::vector<vtkIdType>(numStepsPlusCrest));
    if (numFoldPoints * numStepsPlusCrest != upIds.size()) {
        throw std::runtime_error("Got an unexpected amount of up ids: "
            + std::to_string(numFoldPoints * numStepsPlusCrest) + " != " + std::to_string(upIds.size()));
    }

    double p[3]; // buffer to store points in
    for (size_t l = 0; l < newUpIds.size(); ++l)
    {
        for (size_t s = 0; s < newUpIds[l].size(); ++s)
        {
            spokePoints->GetPoint(upIds[(l * numStepsPlusCrest + s)], p);
            newUpIds[l][s] = spokePoints->InsertNextPoint(p);
            typeArray->InsertNextValue(1);
        }
    }

    // make medial sheet cells
    // for each quad, make two triangles
    //
    //   .-.-.
    //   |/|/|
    //   .-.-.
    for (size_t l = 0; l < newUpIds.size(); ++l)
    {
        auto nextL = (l + 1) % newUpIds.size();
        for (size_t s = 0; s < newUpIds[l].size(); ++s)
        {
            auto nextS = (s + 1) % newUpIds[l].size();
            vtkNew<vtkTriangle> tri1;
            tri1->GetPointIds()->SetId(0, newUpIds[l][s]);
            tri1->GetPointIds()->SetId(1, newUpIds[nextL][s]);
            tri1->GetPointIds()->SetId(2, newUpIds[l][nextS]);

            vtkNew<vtkTriangle> tri2;
            tri2->GetPointIds()->SetId(0, newUpIds[nextL][nextS]);
            tri2->GetPointIds()->SetId(1, newUpIds[nextL][s]);
            tri2->GetPointIds()->SetId(2, newUpIds[l][nextS]);

            medialSheet->InsertNextCell(tri1);
            medialSheet->InsertNextCell(tri2);
        }
    }
}

vtkSmartPointer<vtkPolyData> ReadEllipticalSRepJson(rapidjson::Value& json)
{
    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    vtkNew<vtkPoints> spokePoints;
    vtkNew<vtkCellArray> spokeLines;
    vtkNew<vtkCellArray> medialSheet;
    vtkNew<vtkIntArray> typeArray;
    typeArray->SetName("cellType");
    typeArray->SetNumberOfComponents(1);
    polyData->SetPoints(spokePoints);
    polyData->SetLines(spokeLines);
    polyData->SetPolys(medialSheet);
    polyData->GetPointData()->SetActiveScalars("cellType");
    polyData->GetPointData()->SetScalars(typeArray);

    const auto numFoldPoints = readUint(SafeFindMember(json, keys::CrestPoints)->value);
    const auto numStepsPlusCrest = 1 + readUint(SafeFindMember(json, keys::Steps)->value);
    auto& skeleton = SafeFindMember(json, keys::Skeleton)->value;
    if (!skeleton.IsArray())
    {
        throw std::invalid_argument("Expected a JSON array.");
    }

    std::vector<vtkIdType> upIds;
    std::vector<vtkIdType> crestIds;

    for (auto& row : skeleton.GetArray())
    { //row is the line out from the spine
        if (!row.IsArray())
        {
            throw std::runtime_error("Error parsing in vtkMRMLEllipticalSRepNode. Row is not array.");
        }
        for (auto& object : row.GetArray())
        {
            upIds.push_back(addSpoke(*polyData, SafeFindMember(object, keys::UpSpoke)->value, 0));
            addSpoke(*polyData, SafeFindMember(object, keys::DownSpoke)->value, 4);

            auto crestIter = object.FindMember(keys::CrestSpoke);
            if (crestIter != object.MemberEnd())
            {
                crestIds.push_back(addSpoke(*polyData, crestIter->value, 2));
            }
        }
    }

    addMedialSheet(*polyData, upIds, numFoldPoints, numStepsPlusCrest);
    addCrestCurve(*polyData, crestIds);

    return polyData;
}

} // namespace {}

vtkSmartPointer<vtkPolyData> ShapePopulationData::ReadSRepJson(const std::string& a_filePath)
{
    FileRAII file(a_filePath);
    if (!file.GetHandle())
    {
        return nullptr;
    }

    const size_t bufferSize = 65535;
    std::vector<char> buffer(bufferSize);
    rapidjson::FileReadStream fs(file.GetHandle(), buffer.data(), buffer.size());
    std::unique_ptr<rapidjson::Document> jsonRoot(new rapidjson::Document);
    if (jsonRoot->ParseStream(fs).HasParseError()) {
        return nullptr;
    }

    try
    {
        if (jsonRoot->HasMember(keys::EllipticalSRep))
        {
            return ReadEllipticalSRepJson((*jsonRoot)[keys::EllipticalSRep]);
        }
        return nullptr;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error reading SRep: " << e.what() << std::endl;
        return nullptr;
    }
    catch (...)
    {
        std::cerr << "Unknown error reading SRep" << std::endl;
        return nullptr;
    }
}
