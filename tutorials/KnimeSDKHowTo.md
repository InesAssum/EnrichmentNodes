# Setting up KNIME SDK with Generic KNIME Nodes (GKN)
## Setting up [KNIME SDK](https://github.com/knime/knime-sdk-setup)
1. Download and install [JDK](https://www.oracle.com/technetwork/java/javase/downloads/jdk11-downloads-5066655.html)  and [Eclipse for RCP and RDP Developers](https://www.eclipse.org/downloads/packages/) in case you haven’t already.
2. Clone [KNIME SDK](https://github.com/knime/knime-sdk-setup) into your Eclipse workspace:  
Go to File → Import → Git → Projects from Git File → Clone URI. Enter https://github.com/knime/knime-sdk-setup as URI and proceed.  
Next, select the latest release (releases/2018-12). Next, next, next  and press Finish. 
3. In the imported org.knime.sdk.setup project, you find three target platform definition files (ending with .target).  
A target definition defines the set of KNIME Extensions and Integrations which will be available when starting your KNIME Analytics Platform development version. Double click KNIME-AP-complete.target and click Set as Active Target Platform. Note: Resolving the targets for the first time may take a while (>10 min).

## Importing GenericKnimeNodes to Eclipse Workspace:
1. Clone or download [GKN](https://github.com/genericworkflownodes/GenericKnimeNodes) to your local hard drive in case you haven’t already.
2. From Eclipse Select File -> Open Projects from filesystem. Navigate to your GenericKnimeNodes directory and import com.genericworkflownodes.knime and com.genericworkflownodes.knime.config. 

## Importing Custom Nodes:
1. From Eclipse Select File -> Open Projects from filesystem. Navigate to your GenericKnimeNodes/generated_plugins directory and import the root folder of the latest EnrichmentNodes.  
Note: The Validation may show problems concerning plugins built for different platforms which can be ignored. The Eclipse workspace should now look like this.
![Final Eclipse workspace](http://ascgitlab.helmholtz-muenchen.de/ines.assum/EnrichmentNodes/raw/2332eb559b374f9df0ae1862ff337dc42f45b21e/tutorials/imgs/EclipseWorkspace.png)  
2. Finally, start your KNIME SDK instance by right clicking KNIME-AP-complete.target -> Run as -> Run configurations... :  
Select org.knime.product.KNIME_PRODUCT from the run a product dropdown menu.  
![RunAs](http://ascgitlab.helmholtz-muenchen.de/ines.assum/EnrichmentNodes/raw/2332eb559b374f9df0ae1862ff337dc42f45b21e/tutorials/imgs/RunAsKnime.png)
Click Run. 
 You can now use KNIME SDK. Note that to install new software in KNIME SDK you have to add it to target definition in Eclipse. 