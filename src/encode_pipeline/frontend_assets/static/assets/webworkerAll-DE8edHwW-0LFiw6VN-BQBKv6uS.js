import{w as v,B as G,$ as w,l as B,V as O,a7 as m,p as F,o as P,j as I,d as A,al as C,u as z,a9 as U}from"./GenomeBrowserPanel-DBc1KYog.js";import{b as V}from"./Filter-DNf-FDip-Bj668h4j-WcC7yzHl.js";import"./index-EN7eb3gx.js";var D=`in vec2 aPosition;
out vec2 vTextureCoord;

uniform vec4 uInputSize;
uniform vec4 uOutputFrame;
uniform vec4 uOutputTexture;

vec4 filterVertexPosition( void )
{
    vec2 position = aPosition * uOutputFrame.zw + uOutputFrame.xy;
    
    position.x = position.x * (2.0 / uOutputTexture.x) - 1.0;
    position.y = position.y * (2.0*uOutputTexture.z / uOutputTexture.y) - uOutputTexture.z;

    return vec4(position, 0.0, 1.0);
}

vec2 filterTextureCoord( void )
{
    return aPosition * (uOutputFrame.zw * uInputSize.zw);
}

void main(void)
{
    gl_Position = filterVertexPosition();
    vTextureCoord = filterTextureCoord();
}
`,E=`in vec2 vTextureCoord;
out vec4 finalColor;
uniform sampler2D uTexture;
void main() {
    finalColor = texture(uTexture, vTextureCoord);
}
`,y=`struct GlobalFilterUniforms {
  uInputSize: vec4<f32>,
  uInputPixel: vec4<f32>,
  uInputClamp: vec4<f32>,
  uOutputFrame: vec4<f32>,
  uGlobalFrame: vec4<f32>,
  uOutputTexture: vec4<f32>,
};

@group(0) @binding(0) var <uniform> gfu: GlobalFilterUniforms;
@group(0) @binding(1) var uTexture: texture_2d<f32>;
@group(0) @binding(2) var uSampler: sampler;

struct VSOutput {
  @builtin(position) position: vec4<f32>,
  @location(0) uv: vec2<f32>
};

fn filterVertexPosition(aPosition: vec2<f32>) -> vec4<f32>
{
    var position = aPosition * gfu.uOutputFrame.zw + gfu.uOutputFrame.xy;

    position.x = position.x * (2.0 / gfu.uOutputTexture.x) - 1.0;
    position.y = position.y * (2.0 * gfu.uOutputTexture.z / gfu.uOutputTexture.y) - gfu.uOutputTexture.z;

    return vec4(position, 0.0, 1.0);
}

fn filterTextureCoord(aPosition: vec2<f32>) -> vec2<f32>
{
    return aPosition * (gfu.uOutputFrame.zw * gfu.uInputSize.zw);
}

@vertex
fn mainVertex(
  @location(0) aPosition: vec2<f32>,
) -> VSOutput {
  return VSOutput(
   filterVertexPosition(aPosition),
   filterTextureCoord(aPosition)
  );
}

@fragment
fn mainFragment(
  @location(0) uv: vec2<f32>,
) -> @location(0) vec4<f32> {
    return textureSample(uTexture, uSampler, uv);
}
`;class M extends V{constructor(){const e=C.from({vertex:{source:y,entryPoint:"mainVertex"},fragment:{source:y,entryPoint:"mainFragment"},name:"passthrough-filter"}),t=z.from({vertex:D,fragment:E,name:"passthrough-filter"});super({gpuProgram:e,glProgram:t})}}class S{constructor(e){this._renderer=e}push(e,t,r){this._renderer.renderPipes.batch.break(r),r.add({renderPipeId:"filter",canBundle:!1,action:"pushFilter",container:t,filterEffect:e})}pop(e,t,r){this._renderer.renderPipes.batch.break(r),r.add({renderPipeId:"filter",action:"popFilter",canBundle:!1})}execute(e){e.action==="pushFilter"?this._renderer.filter.push(e):e.action==="popFilter"&&this._renderer.filter.pop()}destroy(){this._renderer=null}}S.extension={type:[v.WebGLPipes,v.WebGPUPipes,v.CanvasPipes],name:"filter"};const k=new P;function W(g,e){e.clear();const t=e.matrix;for(let r=0;r<g.length;r++){const i=g[r];if(i.globalDisplayStatus<7)continue;const n=i.renderGroup??i.parentRenderGroup;n!=null&&n.isCachedAsTexture?e.matrix=k.copyFrom(n.textureOffsetInverseTransform).append(i.worldTransform):n!=null&&n._parentCacheAsTextureRenderGroup?e.matrix=k.copyFrom(n._parentCacheAsTextureRenderGroup.inverseWorldTransform).append(i.groupTransform):e.matrix=i.worldTransform,e.addBounds(i.bounds)}return e.matrix=t,e}const X=new G({attributes:{aPosition:{buffer:new Float32Array([0,0,1,0,1,1,0,1]),format:"float32x2",stride:2*4,offset:0}},indexBuffer:new Uint32Array([0,1,2,0,2,3])});class q{constructor(){this.skip=!1,this.inputTexture=null,this.backTexture=null,this.filters=null,this.bounds=new U,this.container=null,this.blendRequired=!1,this.outputRenderSurface=null,this.globalFrame={x:0,y:0,width:0,height:0},this.firstEnabledIndex=-1,this.lastEnabledIndex=-1}}class R{constructor(e){this._filterStackIndex=0,this._filterStack=[],this._filterGlobalUniforms=new B({uInputSize:{value:new Float32Array(4),type:"vec4<f32>"},uInputPixel:{value:new Float32Array(4),type:"vec4<f32>"},uInputClamp:{value:new Float32Array(4),type:"vec4<f32>"},uOutputFrame:{value:new Float32Array(4),type:"vec4<f32>"},uGlobalFrame:{value:new Float32Array(4),type:"vec4<f32>"},uOutputTexture:{value:new Float32Array(4),type:"vec4<f32>"}}),this._globalFilterBindGroup=new O({}),this.renderer=e}get activeBackTexture(){var e;return(e=this._activeFilterData)==null?void 0:e.backTexture}push(e){const t=this.renderer,r=e.filterEffect.filters,i=this._pushFilterData();i.skip=!1,i.filters=r,i.container=e.container,i.outputRenderSurface=t.renderTarget.renderSurface;const n=t.renderTarget.renderTarget.colorTexture.source,u=n.resolution,o=n.antialias;if(r.every(f=>!f.enabled)){i.skip=!0;return}const d=i.bounds;if(this._calculateFilterArea(e,d),this._calculateFilterBounds(i,t.renderTarget.rootViewPort,o,u,1),i.skip)return;const s=this._getPreviousFilterData(),a=this._findFilterResolution(u);let c=0,l=0;s&&(c=s.bounds.minX,l=s.bounds.minY),this._calculateGlobalFrame(i,c,l,a,n.width,n.height),this._setupFilterTextures(i,d,t,s)}generateFilteredTexture({texture:e,filters:t}){const r=this._pushFilterData();this._activeFilterData=r,r.skip=!1,r.filters=t;const i=e.source,n=i.resolution,u=i.antialias;if(t.every(a=>!a.enabled))return r.skip=!0,e;const o=r.bounds;if(o.addRect(e.frame),this._calculateFilterBounds(r,o.rectangle,u,n,0),r.skip)return e;const d=n;this._calculateGlobalFrame(r,0,0,d,i.width,i.height),r.outputRenderSurface=m.getOptimalTexture(o.width,o.height,r.resolution,r.antialias),r.backTexture=F.EMPTY,r.inputTexture=e,this.renderer.renderTarget.finishRenderPass(),this._applyFiltersToTexture(r,!0);const s=r.outputRenderSurface;return s.source.alphaMode="premultiplied-alpha",s}pop(){const e=this.renderer,t=this._popFilterData();t.skip||(e.globalUniforms.pop(),e.renderTarget.finishRenderPass(),this._activeFilterData=t,this._applyFiltersToTexture(t,!1),t.blendRequired&&m.returnTexture(t.backTexture),m.returnTexture(t.inputTexture))}getBackTexture(e,t,r){const i=e.colorTexture.source._resolution,n=m.getOptimalTexture(t.width,t.height,i,!1);let u=t.minX,o=t.minY;r&&(u-=r.minX,o-=r.minY),u=Math.floor(u*i),o=Math.floor(o*i);const d=Math.ceil(t.width*i),s=Math.ceil(t.height*i);return this.renderer.renderTarget.copyToTexture(e,n,{x:u,y:o},{width:d,height:s},{x:0,y:0}),n}applyFilter(e,t,r,i){const n=this.renderer,u=this._activeFilterData,o=u.outputRenderSurface===r,d=n.renderTarget.rootRenderTarget.colorTexture.source._resolution,s=this._findFilterResolution(d);let a=0,c=0;if(o){const f=this._findPreviousFilterOffset();a=f.x,c=f.y}this._updateFilterUniforms(t,r,u,a,c,s,o,i);const l=e.enabled?e:this._getPassthroughFilter();this._setupBindGroupsAndRender(l,t,n)}calculateSpriteMatrix(e,t){const r=this._activeFilterData,i=e.set(r.inputTexture._source.width,0,0,r.inputTexture._source.height,r.bounds.minX,r.bounds.minY),n=t.worldTransform.copyTo(P.shared),u=t.renderGroup||t.parentRenderGroup;return u&&u.cacheToLocalTransform&&n.prepend(u.cacheToLocalTransform),n.invert(),i.prepend(n),i.scale(1/t.texture.orig.width,1/t.texture.orig.height),i.translate(t.anchor.x,t.anchor.y),i}destroy(){var e;(e=this._passthroughFilter)==null||e.destroy(!0),this._passthroughFilter=null}_getPassthroughFilter(){return this._passthroughFilter??(this._passthroughFilter=new M),this._passthroughFilter}_setupBindGroupsAndRender(e,t,r){if(r.renderPipes.uniformBatch){const i=r.renderPipes.uniformBatch.getUboResource(this._filterGlobalUniforms);this._globalFilterBindGroup.setResource(i,0)}else this._globalFilterBindGroup.setResource(this._filterGlobalUniforms,0);this._globalFilterBindGroup.setResource(t.source,1),this._globalFilterBindGroup.setResource(t.source.style,2),e.groups[0]=this._globalFilterBindGroup,r.encoder.draw({geometry:X,shader:e,state:e._state,topology:"triangle-list"}),r.type===I.WEBGL&&r.renderTarget.finishRenderPass()}_setupFilterTextures(e,t,r,i){if(e.backTexture=F.EMPTY,e.inputTexture=m.getOptimalTexture(t.width,t.height,e.resolution,e.antialias),e.blendRequired){r.renderTarget.finishRenderPass();const n=r.renderTarget.getRenderTarget(e.outputRenderSurface);e.backTexture=this.getBackTexture(n,t,i==null?void 0:i.bounds)}r.renderTarget.bind(e.inputTexture,!0),r.globalUniforms.push({offset:t})}_calculateGlobalFrame(e,t,r,i,n,u){const o=e.globalFrame;o.x=t*i,o.y=r*i,o.width=n*i,o.height=u*i}_updateFilterUniforms(e,t,r,i,n,u,o,d){const s=this._filterGlobalUniforms.uniforms,a=s.uOutputFrame,c=s.uInputSize,l=s.uInputPixel,f=s.uInputClamp,h=s.uGlobalFrame,x=s.uOutputTexture;o?(a[0]=r.bounds.minX-i,a[1]=r.bounds.minY-n):(a[0]=0,a[1]=0),a[2]=e.frame.width,a[3]=e.frame.height,c[0]=e.source.width,c[1]=e.source.height,c[2]=1/c[0],c[3]=1/c[1],l[0]=e.source.pixelWidth,l[1]=e.source.pixelHeight,l[2]=1/l[0],l[3]=1/l[1],f[0]=.5*l[2],f[1]=.5*l[3],f[2]=e.frame.width*c[2]-.5*l[2],f[3]=e.frame.height*c[3]-.5*l[3];const T=this.renderer.renderTarget.rootRenderTarget.colorTexture;h[0]=i*u,h[1]=n*u,h[2]=T.source.width*u,h[3]=T.source.height*u,t instanceof F&&(t.source.resource=null);const b=this.renderer.renderTarget.getRenderTarget(t);this.renderer.renderTarget.bind(t,!!d),t instanceof F?(x[0]=t.frame.width,x[1]=t.frame.height):(x[0]=b.width,x[1]=b.height),x[2]=b.isRoot?-1:1,this._filterGlobalUniforms.update()}_findFilterResolution(e){let t=this._filterStackIndex-1;for(;t>0&&this._filterStack[t].skip;)--t;return t>0&&this._filterStack[t].inputTexture?this._filterStack[t].inputTexture.source._resolution:e}_findPreviousFilterOffset(){let e=0,t=0,r=this._filterStackIndex;for(;r>0;){r--;const i=this._filterStack[r];if(!i.skip){e=i.bounds.minX,t=i.bounds.minY;break}}return{x:e,y:t}}_calculateFilterArea(e,t){if(e.renderables?W(e.renderables,t):e.filterEffect.filterArea?(t.clear(),t.addRect(e.filterEffect.filterArea),t.applyMatrix(e.container.worldTransform)):e.container.getFastGlobalBounds(!0,t),e.container){const r=(e.container.renderGroup||e.container.parentRenderGroup).cacheToLocalTransform;r&&t.applyMatrix(r)}}_applyFiltersToTexture(e,t){const r=e.inputTexture,i=e.bounds,n=e.filters,u=e.firstEnabledIndex,o=e.lastEnabledIndex;if(this._globalFilterBindGroup.setResource(r.source.style,2),this._globalFilterBindGroup.setResource(e.backTexture.source,3),u===o)n[u].apply(this,r,e.outputRenderSurface,t);else{let d=e.inputTexture;const s=m.getOptimalTexture(i.width,i.height,d.source._resolution,!1);let a=s;for(let c=u;c<o;c++){const l=n[c];if(!l.enabled)continue;l.apply(this,d,a,!0);const f=d;d=a,a=f}n[o].apply(this,d,e.outputRenderSurface,t),m.returnTexture(s)}}_calculateFilterBounds(e,t,r,i,n){var u;const o=this.renderer,d=e.bounds,s=e.filters;let a=1/0,c=0,l=!0,f=!1,h=!1,x=!0,T=-1,b=-1;for(let _=0;_<s.length;_++){const p=s[_];if(p.enabled){if(T===-1&&(T=_),b=_,a=Math.min(a,p.resolution==="inherit"?i:p.resolution),c+=p.padding,p.antialias==="off"?l=!1:p.antialias==="inherit"&&l&&(l=r),p.clipToViewport||(x=!1),!(p.compatibleRenderers&o.type)){h=!1;break}if(p.blendRequired&&!(((u=o.backBuffer)==null?void 0:u.useBackBuffer)??!0)){A("Blend filter requires backBuffer on WebGL renderer to be enabled. Set `useBackBuffer: true` in the renderer options."),h=!1;break}h=!0,f||(f=p.blendRequired)}}if(!h){e.skip=!0;return}if(x&&d.fitBounds(0,t.width/i,0,t.height/i),d.scale(a).ceil().scale(1/a).pad((c|0)*n),!d.isPositive){e.skip=!0;return}e.antialias=l,e.resolution=a,e.blendRequired=f,e.firstEnabledIndex=T,e.lastEnabledIndex=b}_popFilterData(){return this._filterStackIndex--,this._filterStack[this._filterStackIndex]}_getPreviousFilterData(){let e,t=this._filterStackIndex-1;for(;t>0&&(t--,e=this._filterStack[t],!!e.skip););return e}_pushFilterData(){let e=this._filterStack[this._filterStackIndex];return e||(e=this._filterStack[this._filterStackIndex]=new q),this._filterStackIndex++,e}}R.extension={type:[v.WebGLSystem,v.WebGPUSystem],name:"filter"};w.add(R);w.add(S);
