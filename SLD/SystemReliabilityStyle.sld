<?xml version="1.0" encoding="UTF-8"?>
<StyledLayerDescriptor xmlns="http://www.opengis.net/sld" xmlns:ogc="http://www.opengis.net/ogc" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" version="1.1.0" xmlns:xlink="http://www.w3.org/1999/xlink" xsi:schemaLocation="http://www.opengis.net/sld http://schemas.opengis.net/sld/1.1.0/StyledLayerDescriptor.xsd" xmlns:se="http://www.opengis.net/se">
  <NamedLayer>
    <se:Name>output_sr</se:Name>
    <UserStyle>
      <se:Name>output_sr</se:Name>
      <se:FeatureTypeStyle>
        <Rule>
          <Name>Probability from 0 to 0.1 ; COLOR: #d5ebba</Name>
          <Filter xmlns="http://www.opengis.net/ogc">
            <And>
              <PropertyIsGreaterThanOrEqualTo>
                <ogc:PropertyName>Prob_Disruption</ogc:PropertyName>
                <Literal>0</Literal>
              </PropertyIsGreaterThanOrEqualTo>
              <PropertyIsLessThan>
                <ogc:PropertyName>Prob_Disruption</ogc:PropertyName>
                <Literal>0.1</Literal>
              </PropertyIsLessThan>
            </And>
          </Filter>
          <se:PolygonSymbolizer>
            <se:Fill>
              <se:SvgParameter name="fill">#d5ebba</se:SvgParameter>
            </se:Fill>
            <se:Stroke>
              <se:SvgParameter name="stroke">#000001</se:SvgParameter>
              <se:SvgParameter name="stroke-width">1</se:SvgParameter>
              <se:SvgParameter name="stroke-linejoin">bevel</se:SvgParameter>
            </se:Stroke>
          </se:PolygonSymbolizer>
        </Rule>
        <Rule>
          <Name>Probability from 0.1 to 0.25 ; COLOR: #fcffa4</Name>
          <Filter xmlns="http://www.opengis.net/ogc">
            <And>
              <PropertyIsGreaterThanOrEqualTo>
                <ogc:PropertyName>Prob_Disruption</ogc:PropertyName>
                <Literal>0.1</Literal>
              </PropertyIsGreaterThanOrEqualTo>
              <PropertyIsLessThan>
                <ogc:PropertyName>Prob_Disruption</ogc:PropertyName>
                <Literal>0.25</Literal>
              </PropertyIsLessThan>
            </And>
          </Filter>
          <se:PolygonSymbolizer>
            <se:Fill>
              <se:SvgParameter name="fill">#fcffa4</se:SvgParameter>
            </se:Fill>
            <se:Stroke>
              <se:SvgParameter name="stroke">#000001</se:SvgParameter>
              <se:SvgParameter name="stroke-width">1</se:SvgParameter>
              <se:SvgParameter name="stroke-linejoin">bevel</se:SvgParameter>
            </se:Stroke>
          </se:PolygonSymbolizer>
        </Rule>
        <Rule>
          <Name>Probability from 0.25 to 0.5 ; COLOR: #f6db4c</Name>
          <Filter xmlns="http://www.opengis.net/ogc">
            <And>
              <PropertyIsGreaterThanOrEqualTo>
                <ogc:PropertyName>Prob_Disruption</ogc:PropertyName>
                <Literal>0.25</Literal>
              </PropertyIsGreaterThanOrEqualTo>
              <PropertyIsLessThan>
                <ogc:PropertyName>Prob_Disruption</ogc:PropertyName>
                <Literal>0.5</Literal>
              </PropertyIsLessThan>
            </And>
          </Filter>
          <se:PolygonSymbolizer>
            <se:Fill>
              <se:SvgParameter name="fill">#f6db4c</se:SvgParameter>
            </se:Fill>
            <se:Stroke>
              <se:SvgParameter name="stroke">#000001</se:SvgParameter>
              <se:SvgParameter name="stroke-width">1</se:SvgParameter>
              <se:SvgParameter name="stroke-linejoin">bevel</se:SvgParameter>
            </se:Stroke>
          </se:PolygonSymbolizer>
        </Rule>
        <Rule>
          <Name>Probability from 0.5 to 0.75 ; COLOR: #fcad11</Name>
          <Filter xmlns="http://www.opengis.net/ogc">
            <And>
              <PropertyIsGreaterThanOrEqualTo>
                <ogc:PropertyName>Prob_Disruption</ogc:PropertyName>
                <Literal>0.5</Literal>
              </PropertyIsGreaterThanOrEqualTo>
              <PropertyIsLessThan>
                <ogc:PropertyName>Prob_Disruption</ogc:PropertyName>
                <Literal>0.75</Literal>
              </PropertyIsLessThan>
            </And>
          </Filter>
          <se:PolygonSymbolizer>
            <se:Fill>
              <se:SvgParameter name="fill">#fcad11</se:SvgParameter>
            </se:Fill>
            <se:Stroke>
              <se:SvgParameter name="stroke">#000001</se:SvgParameter>
              <se:SvgParameter name="stroke-width">1</se:SvgParameter>
              <se:SvgParameter name="stroke-linejoin">bevel</se:SvgParameter>
            </se:Stroke>
          </se:PolygonSymbolizer>
        </Rule>
        <Rule>
          <Name>Probability from 0.75 to 0.8 ; COLOR: #f88311</Name>
          <Filter xmlns="http://www.opengis.net/ogc">
            <And>
              <PropertyIsGreaterThanOrEqualTo>
                <ogc:PropertyName>Prob_Disruption</ogc:PropertyName>
                <Literal>0.75</Literal>
              </PropertyIsGreaterThanOrEqualTo>
              <PropertyIsLessThan>
                <ogc:PropertyName>Prob_Disruption</ogc:PropertyName>
                <Literal>0.8</Literal>
              </PropertyIsLessThan>
            </And>
          </Filter>
          <se:PolygonSymbolizer>
            <se:Fill>
              <se:SvgParameter name="fill">#f88311</se:SvgParameter>
            </se:Fill>
            <se:Stroke>
              <se:SvgParameter name="stroke">#000001</se:SvgParameter>
              <se:SvgParameter name="stroke-width">1</se:SvgParameter>
              <se:SvgParameter name="stroke-linejoin">bevel</se:SvgParameter>
            </se:Stroke>
          </se:PolygonSymbolizer>
        </Rule>
        <Rule>
          <Name>Probability from 0.8 to 0.9 ; COLOR: #cb4149</Name>
          <Filter xmlns="http://www.opengis.net/ogc">
            <And>
              <PropertyIsGreaterThanOrEqualTo>
                <ogc:PropertyName>Prob_Disruption</ogc:PropertyName>
                <Literal>0.8</Literal>
              </PropertyIsGreaterThanOrEqualTo>
              <PropertyIsLessThan>
                <ogc:PropertyName>Prob_Disruption</ogc:PropertyName>
                <Literal>0.9</Literal>
              </PropertyIsLessThan>
            </And>
          </Filter>
          <se:PolygonSymbolizer>
            <se:Fill>
              <se:SvgParameter name="fill">#cb4149</se:SvgParameter>
            </se:Fill>
            <se:Stroke>
              <se:SvgParameter name="stroke">#000001</se:SvgParameter>
              <se:SvgParameter name="stroke-width">1</se:SvgParameter>
              <se:SvgParameter name="stroke-linejoin">bevel</se:SvgParameter>
            </se:Stroke>
          </se:PolygonSymbolizer>
        </Rule>
        <Rule>
          <Name>Probability from 0.9 to 1.0 ; COLOR: #60126e</Name>
          <Filter xmlns="http://www.opengis.net/ogc">
            <And>
              <PropertyIsGreaterThanOrEqualTo>
                <ogc:PropertyName>Prob_Disruption</ogc:PropertyName>
                <Literal>0.9</Literal>
              </PropertyIsGreaterThanOrEqualTo>
              <PropertyIsLessThan>
                <ogc:PropertyName>Prob_Disruption</ogc:PropertyName>
                <Literal>1.0</Literal>
              </PropertyIsLessThan>
            </And>
          </Filter>
          <se:PolygonSymbolizer>
            <se:Fill>
              <se:SvgParameter name="fill">#60126e</se:SvgParameter>
            </se:Fill>
            <se:Stroke>
              <se:SvgParameter name="stroke">#000001</se:SvgParameter>
              <se:SvgParameter name="stroke-width">1</se:SvgParameter>
              <se:SvgParameter name="stroke-linejoin">bevel</se:SvgParameter>
            </se:Stroke>
          </se:PolygonSymbolizer>
        </Rule>
      </se:FeatureTypeStyle>
    </UserStyle>
  </NamedLayer>
</StyledLayerDescriptor>
